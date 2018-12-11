import argparse, sys, gzip, random
import numpy as np

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread

import genomics

from time import sleep


#######################################################################################################################

'''A function that reads from the window queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def freqs_wrapper(windowQueue, resultQueue, genoFormat, sampleData, minData, target, asCounts, keepNanLines = False):
    while True:
        
        windowNumber,window = windowQueue.get() # retrieve window
        
        #make alignment objects
        aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = genoFormat)
        popAlns = dict(zip(sampleData.popNames, [aln.subset(groups=[pop]) for pop in sampleData.popNames]))
        
        #target base at each site
        
        if target == "derived":
            #use last pop as outgroup
            outgroup = sampleData.popNames[-1]
            inAln = aln.subset(groups = sampleData.popNames[:-1])
            baseColumns = np.array([genomics.derivedAllele(inAln.numArray[:,i][inAln.nanMask[:,i]],
                                                           popAlns[outgroup].numArray[:,i][popAlns[outgroup].nanMask[:,i]],
                                                           numeric=True)
                                    for i in xrange(aln.l)]).reshape([aln.l,1])
            
        else:
            #otherwise get minor allele.
                        
            baseColumns = np.array([genomics.minorAllele(aln.numArray[:,i][aln.nanMask[:,i]]) for i in xrange(aln.l)]).reshape([aln.l,1])
        
        goodSites = np.apply_along_axis(lambda(x): ~np.any(np.isnan(x)),1,baseColumns)
        
        #get freqs per pop
        popFreqs = []
        for pop in sampleData.popNames:
            #first find sites with sufficient data
            goodData = popAlns[pop].siteNonNan() >= minData
            sites = np.where(goodSites & goodData)[0]
            baseFreqs = popAlns[pop].siteFreqs(sites, asCounts=asCounts)
            popColumns = baseColumns[sites,:].astype(int)
            popRows = np.repeat(np.arange(len(sites))[:,np.newaxis],popColumns.shape[1], axis = 1)
            targetFreqs =  np.zeros([aln.l, popColumns.shape[1]], dtype=int if asCounts else float)
            if not asCounts: targetFreqs.fill(np.nan)
            if len(sites) >= 1: targetFreqs[sites,:] = baseFreqs[popRows,popColumns]
            popFreqs.append(np.around(targetFreqs, 4))
        
        allFreqs = np.hstack(popFreqs)
        
        if not keepNanLines:
            if not asCounts:
                outSites = np.where(~np.apply_along_axis(np.all, 1, np.isnan(allFreqs)))[0]
            else: outSites = np.where(~np.apply_along_axis(np.all, 1, allFreqs==0))[0]
        else: outSites = range(aln.l)
                
        outArray = np.column_stack(([window.scaffold]*len(outSites),
                                    np.array(window.positions)[outSites].astype(str),
                                    allFreqs[outSites,:].astype(str),))
        
        resultStrings = ["\t".join(row) for row in outArray]
        
        resultQueue.put((windowNumber, resultStrings,))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, verbose):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    while True:
        windowNumber, results = doneQueue.get()
        resultsReceived += 1
        if verbose:
            print >> sys.stderr, "Sorter received window", windowNumber
        if windowNumber == expect:
            writeQueue.put((windowNumber,results))
            if verbose:
                print >> sys.stderr, "Block", windowNumber, "sent to writer"
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,results))
                    if verbose:
                        print >> sys.stderr, "Block", expect, "sent to writer"
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
        if verbose:
            print >> sys.stderr, "\nWriter received window", windowNumber
        for outLine in results:
            out.write(outLine + "\n")
            linesWritten += 1
        resultsWritten += 1

'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        print >> sys.stderr, windowsQueued, "windows queued | ", resultsReceived, "windows analysed | ", resultsWritten, "windows written |", linesWritten, "lines written"

def lineReader(fileObj):
    line = fileObj.readline()
    while len(line) >= 1:
        yield line
        line = fileObj.readline()

#########################################################################################################################


### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-o", "--outFile", help="Output csv file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Data format for output", action = "store", choices = ("phased","diplo","alleles"), default = "phased")

#populations
parser.add_argument("-p", "--population", help="Pop name and optionally sample names (separated by commas)",
                    required = True, action='append', nargs="+", metavar=("popName","[samples]"))
parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)

#Frequency for all bases or a single base
parser.add_argument("--target", help="All or single base frequency (derived assumes last pop is outgroup)", choices = ("minor","derived"),
                    action = "store", default="minor")

parser.add_argument("--asCounts", help="Return drwquencies as counts", action='store_true')

#define ploidy if not 2
parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")

#optional missing data argument
parser.add_argument("--minData", help="Minimum proportion of non-missing data per population", type=float, action = "store", default = 0.1, metavar = "proportion")

#contigs
parser.add_argument("--include", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="File of contigs (one per line)", action='store')

parser.add_argument("--keepNanLines", help="Output lines with no information", action='store_true')


#other
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--verbose", help="Verbose output.", action = "store_true")
parser.add_argument("--windSize", help="Size of windows to process in each thread", type=int, action = "store", default = 10000)
parser.add_argument("--test", help="Test - runs 10 windows", action='store_true')


args = parser.parse_args()

############## parse populations

popNames = []
popInds = []
for p in args.population:
    popNames.append(p[0])
    if len(p) > 1: popInds.append(p[1].split(","))
    else: popInds.append([])

if args.popsFile:
    with open(args.popsFile, "r") as pf: popDict = dict([ln.split() for ln in pf])
    for ind in popDict.keys():
        try: popInds[popNames.index(popDict[ind])].append(ind)
        except: pass

for p in popInds: assert len(p) >= 1, "All populations must be represented by at least one sample."

allInds = list(set([i for p in popInds for i in p]))

if args.ploidy is not None:
    ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allInds)
    assert len(ploidy) == len(allInds), "Incorrect number of ploidy values supplied."
    ploidyDict = dict(zip(allInds,ploidy))
elif args.ploidyFile is not None:
    with open(args.ploidyFile, "r") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
else: ploidyDict = dict(zip(allInds,[2]*len(allInds)))

sampleData = genomics.SampleData(popNames = popNames, popInds = popInds, ploidyDict = ploidyDict)


############################################################################################################################################

#scafs to exclude

if args.exclude:
    with open(args.exclude, "r") as excludeFile:
        scafsToExclude = [line.rstrip() for line in excludeFile]
    print >> sys.stderr, len(scafsToExclude), "scaffolds will be excluded."
else: scafsToExclude = None

if args.include:
    with open(args.include, "r") as includeFile:
        scafsToInclude = [line.rstrip() for line in includeFile]
    print >> sys.stderr, len(scafsToInclude), "scaffolds will be analysed."
else: scafsToInclude = None


#open files

if args.genoFile:
    if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "r")
    else: genoFile = open(args.genoFile, "r")
else: genoFile = sys.stdin


if args.outFile:
    if args.outFile[-3:] == ".gz": outFile = gzip.open(args.outFile, "w")
    else: outFile = open(args.outFile, "w")
else: outFile = sys.stdout

outFile.write("scaffold\tposition\t")
outFile.write("\t".join(popNames) + "\n")

##########################################################################################################################

#counting stat that will let keep track of how far we are
windowsQueued = 0
resultsReceived = 0
resultsWritten = 0
linesWritten = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
windowQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The command should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
for x in range(args.threads):
  worker = Process(target=freqs_wrapper, args = (windowQueue, resultQueue, args.genoFormat, sampleData,
                                                 args.minData, args.target, args.asCounts, args.keepNanLines,))
  worker.daemon = True
  worker.start()
  print >> sys.stderr, "started worker", x


'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,args.verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, outFile, args.verbose,))
worker.daemon = True
worker.start()

'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()

########################################################################################################################

#generate windows and queue

windowGenerator = genomics.nonOverlappingSitesWindows(genoFile, args.windSize,names=sampleData.indNames,
                                                      include = scafsToInclude,exclude = scafsToExclude)


if not args.test:
    for window in windowGenerator:
        windowQueue.put((windowsQueued,window))
        windowsQueued += 1
else:
    for window in windowGenerator:
        windowQueue.put((windowsQueued,window))
        windowsQueued += 1
        if windowsQueued == 10: break


############################################################################################################################################

print >> sys.stderr, "\nWriting final results...\n"
while resultsWritten < windowsQueued:
  sleep(1)

sleep(5)

genoFile.close()
outFile.close()

print >> sys.stderr, "\nDone."

sys.exit()

