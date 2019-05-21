#!/usr/bin/env python

import argparse
import sys
import gzip
import numpy as np
import itertools

import genomics

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread
from time import sleep



####################################################################################################################################


'''A function that reads from the window queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def stats_wrapper(windowQueue, resultQueue, windType, genoFormat, sampleData, minSites, minPerInd, includeSameWithSame,
                  outFormat, roundTo, outputWindowData, addWindowID=False):
    while True:
        nInd = len(sampleData.indNames)
        windowNumber,window = windowQueue.get() # retrieve window
        if windType == "coordinate" or windType == "predefined":
            scaf,start,end,mid,sites = (window.scaffold, window.limits[0], window.limits[1], window.midPos(),window.seqLen())
        else: scaf,start,end,mid,sites = (window.scaffold, window.firstPos(), window.lastPos(),window.midPos(),window.seqLen())
        if sites >= minSites:
            isGood = True
            #make alignment object
            aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = genoFormat)
            if minPerInd and min(aln.seqNonNan()) < minPerInd: isGood = False
            else:
                pairDistDict = aln.indPairDists(includeSameWithSame=includeSameWithSame)
                distMat = np.zeros([nInd,nInd])
                for i,j in itertools.combinations_with_replacement(range(nInd),2):
                    distMat[i,j] = distMat[j,i] = pairDistDict[sampleData.indNames[i]][sampleData.indNames[j]]
        else: isGood = False
        
        if not isGood:
            distMat = np.empty([nInd,nInd])
            distMat.fill(np.NaN)
        if outFormat == "nexus": distMatString = genomics.makeDistMatNexusString(distMat, names=sampleData.indNames, roundTo=roundTo)
        elif outFormat == "phylip": distMatString = genomics.makeDistMatPhylipString(distMat, names=sampleData.indNames, roundTo=roundTo)
        elif outFormat == "raw": distMatString = genomics.makeDistMatString(distMat, roundTo=roundTo) + "\n"
        result = {"main":distMatString}
        if outputWindowData:
            windowData = [] if not addWindowID else [window.ID]
            windowData += [scaf,start,end,mid,sites]
            windowDataString = "\t".join([str(x) for x in windowData]) + "\n"
            result["windows"] = windowDataString
        resultQueue.put((windowNumber, result, isGood))


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
def writer(writeQueue, outs, writeFailedWindows=False):
    global resultsWritten
    global resultsHandled
    while True:
        resNumber,result,isGood = writeQueue.get()
        if verbose:
            print >> sys.stderr, "Writer received result", resNumber
        if isGood or writeFailedWindows:
            for key in result.keys():
                outs[key].write(result[key])
            resultsWritten += 1
        resultsHandled += 1


'''loop that checks stats'''
def checkStats():
  while True:
    sleep(10)
    print >> sys.stderr, windowsQueued, "windows queued", resultsReceived, "results received", resultsWritten, "results written."


####################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("--windType", action = "store", choices = ("sites","coordinate","predefined","cat"), default = "coordinate",
                    help="Type of windows to makem, or concatenate all sites")

parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-s", "--stepSize", help="Step size for sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-m", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites", default = 1)
parser.add_argument("-Mi", "--minPerInd", help="Minumum good sites per individual", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)
parser.add_argument("--windCoords", help="Window coordinates file (scaffold start end)", required = False)

parser.add_argument("--samples", help="Samples to include for individual analysis", nargs="+", action = "store", metavar = "sample names")
parser.add_argument("--includeSameWithSame", action="store_true",
                    help="Include comparisons of each haplotype to itself.")

parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
parser.add_argument("--haploid", help="Alternatively just name samples that are haploid", nargs="+", action = "store", metavar = "sample names")
parser.add_argument("--inferPloidy", help="Ploidy will be inferred in each window (NOT RECOMMENED)", action = "store_true")

parser.add_argument("-g", "--genoFile", help="Input genotypes file", required = False)
parser.add_argument("-o", "--outFile", help="Results file", required = False)
parser.add_argument("--windowDataOutFile", help="Optional window data file", required = False)

parser.add_argument("-f", "--genoFormat", action='store', choices = ("phased","pairs","haplo","diplo"), required = True,
                    help="Format of genotypes in genotypes file")
parser.add_argument("--outFormat", action = "store", choices = ("raw","phylip","nexus"), default = "phylip",
                    help="Format for distance matrix output")

parser.add_argument("--roundTo", help="Round to N decomal places", type=int, action = "store", default=4)

parser.add_argument("--exclude", help="File of scaffolds to exclude", required = False)
parser.add_argument("--include", help="File of scaffolds to analyse", required = False)

parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
parser.add_argument("--verbose", help="Verbose output", action="store_true")
parser.add_argument("--addWindowID", help="Add window name or number as first column", action="store_true")
parser.add_argument("--writeFailedWindows", help="Write output even for windows with too few sites.", action="store_true")


args = parser.parse_args()

#window parameters
minSites = args.minSites

if args.windType == "coordinate":
    assert args.windSize, "Window size must be provided."
    windSize = args.windSize
    stepSize = args.stepSize
    if not stepSize: stepSize = windSize
    assert not args.overlap, "Overlap does not apply to coordinate windows. Use --stepSize instead."
    assert not args.maxDist, "Maximum distance only applies to sites windows."

elif args.windType == "sites":
    assert args.windSize, "Window size (number of sites) must be provided."
    windSize = args.windSize
    overlap = args.overlap
    if not overlap: overlap = 0
    maxDist = args.maxDist
    if not maxDist: maxDist = np.inf
    assert not args.stepSize, "Step size only applies to coordinate windows. Use --overlap instead."
elif args.windType == "predefined":
    assert args.windCoords, "Please provide a file of window coordinates."
    assert not args.overlap, "Overlap does not apply for predefined windows."
    assert not args.maxDist, "Maximum does not apply for predefined windows."
    assert not args.stepSize,"Step size does not apply for predefined windows."
    assert not args.include,"You cannot only include specific scaffolds if using predefined windows."
    assert not args.exclude,"You cannot exclude specific scaffolds if using predefined windows."
    with open(args.windCoords,"r") as wc: windCoords = tuple([(x,int(y),int(z),) for x,y,z in [line.split()[:3] for line in wc]])
else:
    minSites = 1

if not minSites: minSites = windSize

#other
verbose = args.verbose


############## parse samples and populations
allInds = []

if args.samples is not None:
    allInds = list(set(allInds + args.samples))

if len(allInds) == 0:
    with gzip.open(args.genoFile, "r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r") as gf:
        allInds = gf.readline().split()[2:]

if args.ploidy is not None:
    ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allInds)
    assert len(ploidy) == len(allInds), "Incorrect number of ploidy values supplied."
    ploidyDict = dict(zip(allInds,ploidy))
elif args.ploidyFile is not None:
    with open(args.ploidyFile, "r") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
elif args.inferPloidy:
    ploidyDict = dict(zip(allInds,[None]*len(allInds)))
else:
    if args.genoFormat == "haplo": ploidyDict = dict(zip(allInds,[1]*len(allInds)))
    else: ploidyDict = dict(zip(allInds,[2]*len(allInds)))
    if args.haploid:
        for sample in args.haploid: ploidyDict[sample] = 1

sampleData = genomics.SampleData(indNames = allInds, ploidyDict = ploidyDict)

############################################################################################################################################

#open files

if args.genoFile != None: genoFile = gzip.open(args.genoFile, "r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin

outs = {}

if args.outFile: outs["main"] = gzip.open(args.outFile, "w") if args.outFile.endswith(".gz") else open(args.outFile, "w")
else: outs["main"] = sys.stdout

if args.windowDataOutFile:
    outs["windows"] = gzip.open(args.windowDataOutFile, "w") if args.windowDataOutFile.endswith(".gz") else open(args.windowDataOutFile, "w")
    if not args.addWindowID: outs["windows"].write("scaffold,start,end,mid,sites,")
    else: outs["windows"].write("windowID,scaffold,start,end,mid,sites,")
    outputWindowData = True
else:
    outputWindowData = False

############################################################################################################################################

#scafs to exclude (only works for window methods)

if args.exclude:
    scafsFile = open(args.exclude, "rU")
    scafsToExclude = [line.rstrip() for line in scafsFile.readlines()]
    print >> sys.stderr, len(scafsToExclude), "scaffolds will be excluded."
    scafsFile.close()
else:
    scafsToExclude = None

if args.include:
    scafsFile = open(args.include, "rU")
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
for x in range(args.threads):
    worker = Process(target=stats_wrapper, args = (windowQueue, resultQueue, args.windType, args.genoFormat, sampleData, minSites, args.minPerInd,
                                                    args.includeSameWithSame, args.outFormat, args.roundTo, outputWindowData,args.addWindowID))
    worker.daemon = True
    worker.start()
    print >> sys.stderr, "started worker", x


'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, outs, args.writeFailedWindows,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()




##########################################################

if args.windType == "cat":
    window = genomics.parseGenoFile(genoFile, names=sampleData.indNames)
    windowQueue.put((windowsQueued,window))
    windowsQueued += 1
    
else:
    #get windows and analyse
    if args.windType == "coordinate": windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize,
                                                                                sampleData.indNames,
                                                                                include = scafsToInclude,
                                                                                exclude = scafsToExclude)
    elif args.windType == "sites": windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap,
                                                                            maxDist, minSites, sampleData.indNames,
                                                                            include = scafsToInclude,
                                                                            exclude = scafsToExclude)
    else: windowGenerator = genomics.predefinedCoordWindows(genoFile, windCoords, sampleData.indNames)
    
    for window in windowGenerator:
        windowQueue.put((windowsQueued,window))
        windowsQueued += 1

############################################################################################################################################

while resultsHandled < windowsQueued:
  sleep(1)

sleep(5)

genoFile.close()

for out in outs.values(): out.close()

print >> sys.stderr, str(windowsQueued), "windows were tested.\n"
print >> sys.stderr, str(resultsWritten), "results were written.\n"

print >> sys.stderr, "\nDone."

sys.exit()



