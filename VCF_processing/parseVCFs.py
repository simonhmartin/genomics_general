import argparse, gzip, sys
import parseVCF

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread

from time import sleep

##################################################

def parseAndMerge(fileNames, headData, scaffold, start, end, gtFilters, method, skipIndels, missing, ploidy, outSep):
    n = len(fileNames)
    sitesGenerators = [parseVCF.tabixSites(fileNames[x], scaffold, start, end, headData[x].mainHead) for x in range(n)]

    currentSites = []
    for g in sitesGenerators:
        try: currentSites.append(g.next())
        except: currentSites.append(None) 
    
    outLines = []
    
    for pos in xrange(start,end+1):
        filesRepresented = 0
        outObjects = [scaffold, str(pos)]
        for x in xrange(n):
            if currentSites[x] and currentSites[x].POS == pos:
                #get genotypes and add to output
                if not skipIndels or currentSites[x].getType() is not "indel":
                    genotypes = currentSites[x].getGenotypes(gtFilters,asList=True,withPhase=True,missing=missing,allowOnly="ACGT")
                    filesRepresented += 1
                else: genotypes = ["/".join([missing]*ploidy)]*(len(headData[x].sampleNames))
                try: currentSites[x] = sitesGenerators[x].next()
                except: currentSites[x] = None
            else:
                #if not a match, add Ns for this file, and dont read next line
                genotypes = ["/".join([missing]*ploidy)]*(len(headData[x].sampleNames))
            outObjects += genotypes
        #so now we've created the output, but need to decide if we can write it
        if method == "all" or (method == "union" and filesRepresented >= 1) or (method == "intersect" and filesRepresented == n):
            outLines.append(outSep.join(outObjects) + "\n")
    return outLines #and thats it. Move on to the next site in the genome


def parseAndMergeWrapper(inQueue, outQueue, fileNames, headData, gtFilters, method, skipIndels, outSep):
    while True:
        windowNumber,scaffold,start,end = inQueue.get() # retrieve window data
        parsedLines = parseAndMerge(fileNames, headData, scaffold, start, end, gtFilters, method, skipIndels, outSep)
        outQueue.put((windowNumber, parsedLines,))


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
                print >> sys.stderr, "Window", windowNumber, "sent to writer"
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,results))
                    if verbose:
                        print >> sys.stderr, "Window", expect, "sent to writer"
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
            out.write(outLine)
            linesWritten += 1
        resultsWritten += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        print >> sys.stderr, windowsQueued, "windows queued | ", resultsReceived, "windows parsed | ", resultsWritten, "windows written |", linesWritten, "lines written"

################################################################################################################
### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Input vcf file", action = "append", required = True)
parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")
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

args = parser.parse_args()

infiles = args.infile
outfile = args.outfile

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
    print >> sys.stderr, len(include), "contigs will be included."
    
if len(exclude) >= 1:
    exclude = set(exclude)
    print >> sys.stderr, len(exclude), "contigs will be excluded."

gtFilters = []
if args.gtf:
    for gtf in args.gtf:
        try:
            gtfDict = dict([tuple(i.split("=")) for i in gtf])
            for key in gtfDict.keys():
                assert key in ["flag","min","max", "siteTypes", "gtTypes", "samples"]
            for x in ["siteTypes", "gtTypes", "samples"]:
                if x in gtfDict.keys(): gtfDict[x] = gtfDict[x].split(",")
                else: gtfDict[x] = None
            for x in ["min", "max"]:
                if x not in gtfDict.keys(): gtfDict[x] = None
            gtFilters.append(parseVCF.gtFilter(gtfDict["flag"], gtfDict["min"], gtfDict["max"], gtfDict["samples"], gtfDict["siteTypes"], gtfDict["gtTypes"]))
        except:
            print >> sys.stderr, "Bad genotype filter specification. See help."  
            raise


outSep = " " if args.space else "\t"

verbose = args.verbose
##########################################################################################################################

###########################################################################################################################
### open files

headData = [parseVCF.getHeadData(f) for f in infiles]
samples = [s for ss in [h.sampleNames for h in headData] for s in ss]

if outfile: out = gzip.open(outfile, "w") if outfile.endswith(".gz") else open(outfile, "w")
else: out = sys.stdout

###parse fai file

if args.fai:
    with open(args.fai, "r") as fai: scafLens = [(s,int(l)) for s,l in [ln.split()[:2] for ln in fai]]
    scafs = [x[0] for x in scafLens]
    scafLens = dict(scafLens)

if not args.fai:
    scafs = headData[0].contigs
    scafLens = headData[0].contigLengths

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
for x in range(args.threads):
    worker = Process(target=parseAndMergeWrapper,args=(inQueue, outQueue, infiles, headData, gtFilters, args.method,
                                                       args.skipIndels, args.missing, args.ploidy, outSep))
    worker.daemon = True
    worker.start()

'''start two threads for sorting and writing the results'''
worker = Thread(target=sorter, args=(outQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start one Process for sorting and writing the results'''
worker = Thread(target=writer, args=(writeQueue,out,verbose,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()

##########################################################################################################################

out.write(outSep.join(["#CHROM", "POS"] + samples) + "\n")

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


print >> sys.stderr, "\nWriting final results...\n"
while resultsWritten < windowsQueued:
  sleep(1)

sleep(5)

out.close()

print >> sys.stderr, "\nDone."

out.close()
