import argparse
import sys
import gzip
import numpy as np

import genomics

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread
from time import sleep



####################################################################################################################################


'''A function that reads from the window queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def ABBABABA_wrapper(windowQueue, resultQueue, windType, genoFormat, sampleData, P1, P2, P3, O, minData, minSites, stats = ["ABBA","BABA","D","fd"]):
  while True:
    windowNumber,window = windowQueue.get() # retrieve window
    if windType == "coordinate" or windType == "predefined":
        scaf,start,end,mid,sites = (window.scaffold, window.start, window.end, window.midPos(),window.seqLen())
    else: scaf,start,end,mid,sites = (window.scaffold, window.firstPos(), window.lastPos(),window.midPos(),window.seqLen())
    sitesUsed = np.NaN
    if sites >= minSites:
        #make alignment object
        Aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = genoFormat)
        statsDict = genomics.ABBABABA(Aln, P1, P2, P3, O, minData)
        sitesUsed = statsDict["sitesUsed"]
        if sitesUsed >= minSites:
            isGood = True
            values = [round(statsDict[stat],4) for stat in stats]
        else:
            isGood = False
            values = [np.NaN]*len(stats)
    else:
        isGood = False
        values = [np.NaN]*len(stats)
    resultString = ",".join([str(x) for x in [scaf,start,end,mid,sites,sitesUsed] + values])
    resultQueue.put((windowNumber, resultString, isGood))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the result number, and this increases consecutively.'''
def sorter(resultQueue, writeQueue, verbose):
  global resultsReceived
  sortBuffer = {}
  expect = 0
  while True:
    resNumber,result,isGood = resultQueue.get()
    resultsReceived += 1
    if verbose:
      print >> sys.stderr, "Sorter received result", resNumber
    if resNumber == expect:
      writeQueue.put((resNumber,result,isGood))
      if verbose:
        print >> sys.stderr, "Result", resNumber, "sent to writer"
      expect +=1
      #now check buffer for further results
      while True:
        try:
          result,isGood = sortBuffer.pop(str(expect))
          writeQueue.put((expect,result,isGood))
          if verbose:
            print >> sys.stderr, "Result", expect, "sent to writer"
          expect +=1
        except:
          break
    else:
      #otherwise this line is ahead of us, so add to buffer dictionary
      sortBuffer[str(resNumber)] = (result,isGood)

'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out):
  global resultsWritten
  global resultsHandled
  while True:
    resNumber,result,isGood = writeQueue.get()
    if verbose:
      print >> sys.stderr, "Writer received result", resNumber
      if isGood:
        print >> sys.stderr, "Writing good result."
      else:
        print >> sys.stderr, "Omitting bad result."
    if isGood:
      out.write(result + "\n")
      resultsWritten += 1
    resultsHandled += 1


'''loop that checks stats'''
def checkStats():
  while True:
    sleep(10)
    print >> sys.stderr, windowsQueued, "windows queued", resultsReceived, "results received", resultsWritten, "results written."


####################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("--windType", help="Type of windows to make", action = "store", choices = ("sites","coordinate","predefined"), default = "coordinate")
parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-s", "--stepSize", help="Step size for sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-m", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites", default = 1)
parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)
parser.add_argument("--windCoords", help="Window coordinates file (scaffold start end)", required = False)
parser.add_argument("--minData", help="Min proportion of samples genotped per site", type=float, action="store", required = False, default = 0.01, metavar = "proportion")

parser.add_argument("-p", "--population", help="Pop name and ind names (separated by commas)", action='append', nargs=2, required = True, metavar=("name","inds"))
parser.add_argument("--haploid", help="Samples that are haploid (comma separated)", action = "store", metavar = "sample names")

parser.add_argument("-g", "--genoFile", help="Input genotypes file", required = True)
parser.add_argument("-o", "--outFile", help="Results file", required = True)
parser.add_argument("--exclude", help="File of scaffolds to exclude", required = False)
parser.add_argument("--include", help="File of scaffolds to analyse", required = False)
parser.add_argument("-f", "--genoFormat", help="Format of genotypes in genotypes file", action='store', choices = ("phased","pairs","haplo","diplo"), required = True)

parser.add_argument("-T", "--Threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
parser.add_argument("--verbose", help="Verbose output", action="store_true")


args = parser.parse_args()

#window parameters
windType = args.windType

if args.windType == "coordinate":
    assert args.windSize, "Window size must be provided."
    windSize = args.windSize
    stepSize = args.stepSize
    if not stepSize: stepSize = windSize
    assert not args.overlap, "Overlap does not apply to coordinate windows. Use --stepSize instead."
    assert not args.maxDist, "Maximum distance only applies to sites windows."

elif args.windType == "sites":
    assert args.windSize, "Window size (number of sites) must be provided."
    overlap = args.overlap
    if not overlap: overlap = 0
    maxDist = args.maxDist
    if not maxDist: maxDist = np.inf
    assert not args.stepSize, "Step size only applies to coordinate windows. Use --overlap instead."
else:
    assert args.windCoords, "Please provide a file of window coordinates."
    assert not args.overlap, "Overlap does not apply for predefined windows."
    assert not args.maxDist, "Maximum does not apply for predefined windows."
    assert not args.stepSize,"Step size does not apply for predefined windows."
    assert not args.include,"You cannot only include specific scaffolds if using predefined windows."
    assert not args.exclude,"You cannot exclude specific scaffolds if using predefined windows."
    with open(args.windCoords,"r") as wc: windCoords = tuple([(x,int(y),int(z),) for x,y,z in [line.split()[:3] for line in wc]])

minSites = args.minSites
if not minSites: minSites = windSize


minData = args.minData
assert 0 <= minData <= 1, "minimum data per site must be between 0 and 1."

#file info
genoFileName = args.genoFile
genoFormat = args.genoFormat

outFileName = args.outFile

exclude = args.exclude
include = args.include

#other
threads = args.Threads
verbose = args.verbose


############## parse populations
pops = args.population

popNames = [p[0] for p in pops]
print popNames
assert sorted(popNames) == ["O","P1","P2","P3"], "Populations must be named P1, P2, P3 and O."

popInds = [p[1].split(",") for p in pops]

allInds = list(set([i for p in popInds for i in p]))

ploidyDict = dict(zip(allInds,[2]*len(allInds)))
if args.haploid:
    for sample in args.haploid.split(","):
        ploidyDict[sample] = 1

sampleData = genomics.SampleData(popNames = popNames, popInds = popInds, ploidyDict = ploidyDict)

############################################################################################################################################

#open files

if genoFileName[-3:] == ".gz":
  genoFile = gzip.open(genoFileName, "r")
else:
  genoFile = open(genoFileName, "r")

outFile = open(outFileName, "w")

outFile.write("scaffold,start,end,mid,sites,sitesUsed,ABBA,BABA,D,fd\n")

##############################################################

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
  worker = Process(target=ABBABABA_wrapper, args = (windowQueue, resultQueue, windType, genoFormat, sampleData, "P1", "P2", "P3", "O", minData, minSites,))
  worker.daemon = True
  worker.start()


'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, outFile,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()




##########################################################

#get windows and analyse
if windType == "coordinate": windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize,
                                                                            sampleData.indNames,
                                                                            include = scafsToInclude,
                                                                            exclude = scafsToExclude)
elif windType == "sites": windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap,
                                                                         maxDist, minSites, sampleData.indNames,
                                                                         include = scafsToInclude,
                                                                         exclude = scafsToExclude)
else: windowGenerator = genomics.predefinedCoordWindows(genoFile, windCoords, sampleData.indNames)


for window in windowGenerator:
    windowQueue.put((windowsQueued,window))
    windowsQueued += 1

############################################################################################################################################

print >> sys.stderr, "\nWriting final results...\n"
while resultsHandled < windowsQueued:
  sleep(1)

sleep(5)

genoFile.close()
outFile.close()

print >> sys.stderr, str(windowsQueued), "windows were tested.\n"
print >> sys.stderr, str(resultsWritten), "results were written.\n"

print "\nDone."

sys.exit()



