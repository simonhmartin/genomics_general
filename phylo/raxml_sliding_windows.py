
import argparse, sys, os, gzip
import numpy as np
import genomics

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread
from time import sleep



####################################################################################################################################

#command that actually runs raxml and gets the bestTree 
def raxTree(seqBlock, indNames, model, raxml, outgroup = None, uniqueTag = "", test = False, log="/dev/null"):
    #write file
    tempAln = open("temp." + uniqueTag + ".phy", "w")
    tempAln.write(genomics.makeAlnString(indNames,seqBlock))
    tempAln.close()
    if outgroup is not None:
        og = " -o " + ",".join(outgroup)
    else:
        og = ""
    #raxCommand = raxml + " -s temp." + uniqueTag + ".phy -n " + uniqueTag + " -m " + model + og + " -V -f d -p 12345 --silent"
    raxCommand = raxml + " -s temp." + uniqueTag + ".phy -n " + uniqueTag + " -m " + model + og + " -V -f d -p 12345 --silent >>" + log
    if test: print >> sys.stderr, "raxml command:\n", raxCommand
    os.system(raxCommand)
    tempAln.close()
    #try retrieve the result  
    try:
        treeFile = open("RAxML_bestTree." + uniqueTag, "r")
        tree = treeFile.readline()
        treeFile.close()
    except:
        tree = "NA\n"
    #remove files
    if not test:
        os.system("rm temp." + uniqueTag + ".phy*")
        os.system("rm RAxML*" + uniqueTag)
    return tree


'''A function that reads from the window queue, calls sume other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using'''
def raxml_wrapper(windowQueue, resultQueue, windType, genoFormat, model, outgroup, raxml, minSites, minPerInd, test = False):
    while True:
        windowNumber,window = windowQueue.get()
        sites = window.seqLen()
        if test or verbose:
            print >> sys.stderr, "Window", windowNumber, "received for analysis, length:", sites
        if windType == "coordinate": scaf,start,end,mid = (window.scaffold, window.start, window.end, window.midPos())
        else: scaf,start,end,mid = (window.scaffold, window.firstPos(), window.lastPos(), window.midPos())
        data = [window.scaffold, str(start), str(end), str(mid), str(sites)]
        if sites >= minSites:
            block = window.seqs
            indNames = window.names
            #convert block if necessary
            if genoFormat == "phased":
                block = [seq for phasedSeq in [genomics.parsePhase(x) for x in block] for seq in phasedSeq]
                indNames = [i + "_" + x for i in indNames for x in ("A","B")]
            sitesPerInd = [len([x for x in seq if x != "N"]) for seq in block]
            if min(sitesPerInd) >= minPerInd:
                tree = raxTree(block,indNames,model,raxml, outgroup, uniqueTag = scaf + "_" + str(start), test = test, log = log)
            else: tree= "NA\n"
        else: tree = "NA\n"
        
        resultQueue.put((windowNumber, "\t".join(data),tree))


'''a function that watches the result queue and sorts results.'''
def sorter(resultQueue, writeQueue, verbose):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    while True:
        resNumber,data,tree = resultQueue.get()
        resultsReceived += 1
        if verbose:
            print >> sys.stderr, "Sorter received result", resNumber
        if resNumber == expect:
            writeQueue.put((resNumber,data,tree))
            if verbose:
                print >> sys.stderr, "Result", resNumber, "sent to writer"
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    data,tree = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,data,tree))
                    if verbose:
                        print >> sys.stderr, "Result", expect, "sent to writer"
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(resNumber)] = (data,tree)

'''a writer function that writes the sorted result.'''
def writer(writeQueue, dataFile, treesFile):
    global resultsWritten
    global resultsHandled
    while True:
        resNumber,data,tree = writeQueue.get()
        if verbose:
            print >> sys.stderr, "Writer received result", resNumber
        dataFile.write(data + "\n")
        treesFile.write(tree)
        resultsWritten += 1
        resultsHandled += 1


'''loop that checks stats'''
def checkStats():
    while True:
        sleep(10)
        print >> sys.stderr, windowsQueued, "windows queued | ", resultsReceived, "results received | ", resultsWritten, "results written."


####################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("--windType", help="Type of windows to make", action = "store", choices = ("sites","coordinate"), default = "coordinate")

parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = True, metavar="sites")
parser.add_argument("-M", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-Mi", "--minPerInd", help="Minumum good sites per individual", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-S", "--stepSize", help="Step size for coordinate sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)

parser.add_argument("-g", "--genoFile", help="Input genotypes file", required = True)
parser.add_argument("-p", "--prefix", help="Prefix for output files", required = True)

parser.add_argument("--exclude", help="File of scaffolds to exclude", required = False)
parser.add_argument("--include", help="File of scaffolds to analyse", required = False)

parser.add_argument("-f", "--genoFormat", help="Format of genotypes in genotypes file", action='store', choices = ("phased","haplo","diplo"), required = True)

parser.add_argument("--individuals", help="Individuals to include, separated by comma", action = "store", metavar="ind1,ind2,ind3...")

parser.add_argument("--outgroup", help="Outgroup individuals, separated by comma", action = "store", metavar="ind1,ind2,ind3...")

parser.add_argument("--raxml", help="path to raxml executable", action = "store", metavar="path/to/raxml", default="raxml")

parser.add_argument("--model", help="RAxML model", action = "store", default="GTRCAT")

parser.add_argument("--log", help="raxml log file, if you want one.", action = "store", metavar="path/to/log", default="/dev/null")

parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
parser.add_argument("--test", help="run test", action="store_true")
parser.add_argument("--verbose", help="Verbose output", action="store_true")


args = parser.parse_args()

verbose = args.verbose
test = args.test

windType = args.windType

windSize = args.windSize

if args.windType == "coordinate":
    stepSize = args.stepSize
    if not stepSize: stepSize = windSize
    assert not args.overlap, "Overlap noes not apply to coordinate windows. Use --stepSize instead."
    assert not args.maxDist, "Maximum distance only applies to sites windows."

else:
    overlap = args.overlap
    if not overlap: overlap = 0
    maxDist = args.maxDist
    if not maxDist: maxDist = np.inf
    assert not args.stepSize, "Step size only applies to coordinate windows. Use --overlap instead."
    

minSites = args.minSites
if not minSites: minSites = windSize

minPerInd = args.minPerInd
if not minPerInd: minPerInd = minSites


genoFileName = args.genoFile
genoFormat = args.genoFormat


if args.individuals: indNames = args.individuals.split(",")
else: indNames = None

if args.outgroup:
    outgroup = args.outgroup.split(",")
    if genoFormat == "phased":
        outgroup = [i + "_" + x for i in outgroup for x in ("A","B")]
    if test or verbose: print >> sys.stderr, "outgroups:", " ".join(outgroup)

else: outgroup = None


prefix = args.prefix

log = args.log

exclude = args.exclude
include = args.include

threads = args.threads

raxml = args.raxml
model = args.model

############################################################################################################################################

#open files

if args.genoFile: genoFile = gzip.open(args.genoFile, "r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin


dataFile = open(prefix + ".data.tsv", "w")
dataFile.write("scaffold\tstart\tend\tmid\tsites\n")

treesFile = gzip.open(prefix + ".trees.gz", "w") 

############################################################################################################################################


#scafs to exclude

if exclude:
    scafsFile = open(exclude, "rU")
    scafsToExclude = [line.rstrip() for line in scafsFile.readlines()]
    print >> sys.stderr, len(scafsToExclude), "scaffolds will be excluded."
    scafsFile.close()
else:
    scafsToExclude = None

if include:
    scafsFile = open(include, "rU")
    scafsToInclude = [line.rstrip() for line in scafsFile.readlines()]
    print >> sys.stderr, len(scafsToInclude), "scaffolds will be analysed."
    scafsFile.close()
else:
    scafsToInclude = None


##########################################################################################################

#counting stat that will let keep track of how far we are
windowsQueued = 0
resultsReceived = 0
resultsWritten = 0
resultsHandled = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
windowQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''

for x in range(threads):
    worker = Process(target=raxml_wrapper, args = (windowQueue, resultQueue, windType, genoFormat, model, outgroup, raxml, minSites, minPerInd, test,))
    worker.daemon = True
    worker.start()
    print >> sys.stderr, "started worker", x
    

'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, dataFile, treesFile,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()




##########################################################

#get windows and analyse
if windType == "coordinate":
    windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize, indNames, include = scafsToInclude, exclude = scafsToExclude)
else:
    windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap, maxDist, minSites, include = scafsToInclude, exclude = scafsToExclude)
    

for window in windowGenerator:
    #simpleque has no max, so to make sure we haven't gotten ahead of ourselves, we compare windowsQueued to resultsReceived
    while windowsQueued - resultsReceived >= 50:
        sleep(10)
        if test or verbose: print >> sys.stderr, "Waiting for queue to clear..."
    
    if test or verbose:
        print >> sys.stderr, "Sending window", windowsQueued, "to queue. Length:", window.seqLen()
    
    windowQueue.put((windowsQueued,window))
    windowsQueued += 1
    if test and windowsQueued == 10: break

############################################################################################################################################

print >> sys.stderr, "\nWriting final results...\n"
while resultsHandled < windowsQueued:
  sleep(1)

sleep(5)

dataFile.close()
treesFile.close()

print >> sys.stderr, str(windowsQueued), "windows were analysed.\n"
print >> sys.stderr, str(resultsWritten), "results were written.\n"

print "\nDone."

sys.exit()

