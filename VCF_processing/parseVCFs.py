#!/usr/bin/env python

import argparse, gzip, sys, subprocess
import numpy as np
import parseVCF

from multiprocessing import Process

if sys.version_info.major < 3:
    from multiprocessing.queues import SimpleQueue
else:
    from multiprocessing import SimpleQueue

from threading import Thread

from time import sleep

##################################################

#def tabixStream(fileName, chrom, start, end):
    #region = chrom+":"+str(start)+"-"+str(end)  
    #return subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1)

#improved tabix stream that fixes the bug where readining is terminated prematurely
def tabixStream(fileName, region = None, chrom = None, start=None, end=None, header=False):
    if not region: region = chrom+":"+str(start)+"-"+str(end)  
    
    if header:
        p = subprocess.Popen(['tabix', '-h', fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    else:
        p = subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    
    return iter(p.communicate()[0].strip().split("\n"))

def parseAndMerge(fileNames, headData, scaffold, start, end, gtFilters, method, skipIndels, missing, ploidy, outSep, verbose):
    n = len(fileNames)
    
    if verbose: sys.stderr.write("Attempting to extract {}:{}-{}.\n".format(scaffold,start,end))
    
    sitesGenerators = [parseVCF.parseVcfSites(tabixStream(fileNames[x], chrom=scaffold, start=start, end=end),
                                              headData[x]["mainHeaders"], excludeDuplicates=True) for x in range(n)]
    
    currentSites = []
    for x in range(n):
        try: currentSites.append(next(sitesGenerators[x]))
        except:
            sys.stderr.write("WARNING empty window: " + fileNames[x] + " " + scaffold + " " + str(start) + "-" + str(end) + "\n")
            currentSites.append(None) 
    
    outLines = []
    
    sitesConsidered = 0
    linesParsed = 0
    
    for pos in range(start,end+1):
        sitesConsidered+=1
        filesRepresented = 0
        outObjects = [scaffold, str(pos)]
        for x in range(n):
            if currentSites[x] and currentSites[x].POS == pos:
                #get genotypes and add to output
                if not skipIndels or currentSites[x].getType() is not "indel":
                    genotypes = currentSites[x].getGenotypes(gtFilters,asList=True,withPhase=True,missing=missing,allowOnly="ACGT",keepPartial=False)
                    filesRepresented += 1
                else: genotypes = ["/".join([missing]*ploidy)]*headData[x]["nSamples"]
                try: currentSites[x] = next(sitesGenerators[x])
                except: currentSites[x] = None
            else:
                #if not a match, add Ns for this file, and dont read next line
                genotypes = ["/".join([missing]*ploidy)]*headData[x]["nSamples"]
            outObjects += genotypes
        #so now we've created the output, but need to decide if we can write it
        if method == "all" or (method == "union" and filesRepresented >= 1) or (method == "intersect" and filesRepresented == n):
            outLines.append(outSep.join(outObjects) + "\n")
            linesParsed+=1
    if verbose: sys.stderr.write(str(sitesConsidered) + " sites considered. " + str(linesParsed) + " lines parsed.\n")
    return outLines #and thats it. Move on to the next site in the genome


def parseAndMergeWrapper(inQueue, outQueue, fileNames, headData, gtFilters, method, skipIndels, missing, ploidy, outSep, verbose):
    while True:
        
        windowNumber,scaffold,start,end = inQueue.get() # retrieve window data
        
        if windowNumber == -1:
            outQueue.put((-1,None,)) # this is the way of telling everything we're done
            break
        
        parsedLines = parseAndMerge(fileNames, headData, scaffold, start, end, gtFilters, method, skipIndels, missing, ploidy, outSep, verbose)
        outQueue.put((windowNumber, parsedLines,))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, nWorkerThreads, verbose):
    global resultsReceived
    threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
    sortBuffer = {}
    expect = 0
    while True:
        windowNumber, results = doneQueue.get()
        #check if we're finished
        if windowNumber == -1: threadsComplete += 1
        #once all threads report they are done, tel the writer we're done
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1,None,))
            break #this is the way of telling everything we're done
        
        resultsReceived += 1
        if windowNumber == expect:
            writeQueue.put((windowNumber,results))
            if verbose:
                sys.stderr.write("Window {} sent to writer\n".format(windowNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,results))
                    if verbose:
                        sys.stderr.write("Window {} sent to writer\n".format(expect))
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(windowNumber)] = results


'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out, verbose):
    global resultsWritten
    global linesWritten
    while True:
        windowNumber, results = writeQueue.get()
        if windowNumber == -1: break # this is the signal from the sorter that we're done
        if verbose:
            sys.stderr.write("\nWriting window {}\n".format(windowNumber))
        for outLine in results:
            out.write(outLine)
            linesWritten += 1
        resultsWritten += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} windows queued | {} windows parsed | {} windows written | {} lines written\n".format(windowsQueued,
                                                                                                                resultsReceived,
                                                                                                                resultsWritten,
                                                                                                                linesWritten))

################################################################################################################
### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inFile", help="Input vcf file", action = "append", required = True)
parser.add_argument("-o", "--outFile", help="Output csv file", action = "store")
parser.add_argument("-f", "--fai", help="Fasta index file", action = "store")
parser.add_argument("-M", "--method", help="How to merge", action = "store", choices = ("all","intersect","union"), default = "union")

#contigs
parser.add_argument("--include", help="include contigs (separated by commas)", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs (separated by commas)", action='store')
parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')

parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--verbose", help="Verbose output.", action = "store_true")
parser.add_argument("--windSize", help="Size of windows to process in each thread", type=int, action = "store", default = 100000)
parser.add_argument("--space", help="Output separator is space instead of tab", action = "store_true")
parser.add_argument("--test", help="Test - runs 10 windows", action='store_true')


#vcf parsing arguments
parser.add_argument("--gtf", help="Genotype filter. Syntax: flag=X min=X max=X siteTypes=X,X.. gtTypes=X,X.. samples=X,X..", action = "append", nargs = '+')
parser.add_argument("--skipIndels", help="Skip indels", action = "store_true")
parser.add_argument("--missing", help="Value to use for missing data", action = "store", default = "N")
parser.add_argument("--ploidy", help="Ploidy for missing data", action = "store", type = int, default = 2)
parser.add_argument("--outSep", help="Output separator", action = "store", default = "\t")

args = parser.parse_args()

include = []
exclude = []

if args.include: include += args.include.split(",")
if args.exclude: exclude += args.exclude.split(",")

if args.includeFile:
    with open(args.includeFile, 'r') as includeFile:
        include += [c.strip() for c in includeFile.read().split("\n")]

if args.excludeFile:
    with open(args.excludeFile, 'r') as excludeFile:
        exclude += [c.strip() for c in excludeFile.read().split("\n")]

if len(include) >= 1:
    include = set(include)
    sys.stderr.write("{} contigs will be included.".format(len(include)))

if len(exclude) >= 1:
    exclude = set(exclude)
    sys.stderr.write("{} contigs will be excluded.".format(len(exclude)))

gtFilters = [parseVCF.parseGenotypeFilterArg(gtf) for gtf in args.gtf] if args.gtf else []

verbose = args.verbose
##########################################################################################################################

###########################################################################################################################
### open files

headData = [parseVCF.getHeadData(f) for f in args.inFile]
samples = [s for ss in [h["sampleNames"] for h in headData] for s in ss]

if args.outFile: outFile = gzip.open(args.outFile, "w") if args.outFile.endswith(".gz") else open(args.outFile, "w")
else: outFile = sys.stdout

###parse fai file

if args.fai:
    with open(args.fai, "r") as fai: scafLens = [(s,int(l)) for s,l in [ln.split()[:2] for ln in fai]]
    scafs = [x[0] for x in scafLens]
    scafLens = dict(scafLens)

if not args.fai:
    scafs = headData[0]["contigs"]
    scafLens = headData[0]["contigLengths"]

##########################################################################################################

#counting stat that will let keep track of how far we are
windowsQueued = 0
resultsReceived = 0
resultsWritten = 0
linesWritten = 0


'''Create queues to hold the data. One will hold the pod info to be passed to the parser'''
inQueue = SimpleQueue()
#one will hold the results (in the order they come)
outQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for parser. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s)
This one reads from the pod queue, passes each line some analysis function(s), gets the results and sends to the result queue'''
workerThreads = []
sys.stderr.write("\nStarting {} worker threads\n".format(args.threads))
for x in range(args.threads):
    workerThread = Process(target=parseAndMergeWrapper,args=(inQueue, outQueue, args.inFile, headData, gtFilters, args.method,
                                                       args.skipIndels, args.missing, args.ploidy, args.outSep, verbose,))
    workerThread.daemon = True
    workerThread.start()
    workerThreads.append(workerThread)

'''start two threads for sorting and writing the results'''
sorterThread = Thread(target=sorter, args=(outQueue, writeQueue, args.threads, verbose,))
sorterThread.daemon = True
sorterThread.start()

'''start one Process for sorting and writing the results'''
writerThread = Thread(target=writer, args=(writeQueue, outFile, verbose,))
writerThread.daemon = True
writerThread.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
checkerThread = Thread(target=checkStats)
checkerThread.daemon = True
checkerThread.start()

##########################################################################################################################

outFile.write(args.outSep.join(["#CHROM", "POS"] + samples) + "\n")

'''now we go through assuming all files are ordered as in the fai.
if we don't find the line we're looking for we move on to the next'''

for scaf in scafs:
    if (exclude and scaf in exclude) or (include and scaf not in include): continue
    starts = range(1,scafLens[scaf] + 1, args.windSize)
    ends = [s + args.windSize - 1 for s in starts]
    for x in range(len(starts)):
        inQueue.put((windowsQueued,scaf,starts[x],ends[x],))
        windowsQueued += 1
        if args.test and windowsQueued == 10: break
    if args.test and windowsQueued == 10: break


############################################################################################################################################

#Now we send completion signals to all worker threads
for x in range(args.threads):
    inQueue.put((-1,None,None,None,)) # -1 tells the threads to break

#and wait for all to finish
for x in range(len(workerThreads)):
    workerThreads[x].join()

sorterThread.join()
writerThread.join()

sys.stderr.write("\nDone\n")

outFile.close()

sys.exit()
