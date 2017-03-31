#!/usr/bin/env python

import argparse, sys, gzip, random

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread

import genomics

from time import sleep


#######################################################################################################################

'''main worker function. This will watch the inQueue for pods, and pass lines from these pods to be parsed and filtered, before packaging back into a pod and sending on to the resultQueue'''
def analysisWrapper(inQueue,outQueue,inputGenoFormat,outputGenoFormat,headers,include,exclude,samples,minCalls,minPopCalls,
                    minAlleles,maxAlleles,minVarCount,maxHet,minFreq,maxFreq,
                    HWE_P,HWE_side,popDict,ploidyDict,fixed,skipChecks,thinDist):
    while True:
        podNumber,inPod = inQueue.get()
        if verbose: print >> sys.stderr, "Pod", podNumber, "received for analysis."
        outPod = []
        lastScaf = None
        for lineData in inPod:
            lineNumber,line = lineData
            #if verbose: print >> sys.stderr, "Analysing line", lineNumber
            objects = line.split()
            if (include and objects[0] not in include) or (exclude and objects[0] in exclude): continue
            site = genomics.GenomeSite(genotypes=objects[2:], sampleNames=headers[2:], popDict=popDict,
                                       ploidyDict=ploidyDict, genoFormat=inputGenoFormat, skipChecks=skipChecks)
            goodSite = True
            if thinDist:
                pos = int(objects[1])
                if lastScaf != objects[0]:
                    lastPos = pos
                    lastScaf = objects[0]
                    goodSite = False
                elif pos - lastPos < thinDist: goodSite = False
            if goodSite: goodSite = genomics.siteTest(site,samples=samples,minCalls=minCalls,minPopCalls=minPopCalls,
                                               minAlleles=minAlleles,maxAlleles=maxAlleles,minVarCount=minVarCount,
                                               maxHet=maxHet,minFreq=minFreq,maxFreq=maxFreq,HWE_P=HWE_P,HWE_side=HWE_side,fixed=fixed)
            if goodSite:
                outLine = "\t".join(objects[:2] + [str(g) for g in site.asList(samples, mode=outputGenoFormat)]) + "\n"
                outPod.append((lineNumber,outLine))
                if thinDist: lastPos = int(objects[1])
            #if verbose: print >> sys.stderr, objects[0], objects[1], "passed: ", goodSite
        outQueue.put((podNumber,outPod))
        if verbose: print >> sys.stderr, "Pod", podNumber, "analysed, sent to sorter."



'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, verbose):
    global podsDone
    global podsSorted
    sortBuffer = {}
    expect = 0
    while True:
        podNumber, donePod = doneQueue.get()
        podsDone += 1
        if verbose:
            print >> sys.stderr, "Sorter received pod", podNumber
        if podNumber == expect:
            writeQueue.put((podNumber,donePod))
            podsSorted += 1
            if verbose:
                print >> sys.stderr, "Pod", podNumber, "sent to writer"
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    donePod = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,donePod))
                    podsSorted += 1
                    if verbose:
                        print >> sys.stderr, "Pod", expect, "sent to writer"
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(podNumber)] = donePod



'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out):
    global podsWritten
    global linesWritten
    while True:
        podNumber, donePod = writeQueue.get()
        if verbose:
            print >> sys.stderr, "\nWriter received pod", podNumber
        for thing in donePod:
            lineNumber,outLine = thing
            out.write(outLine)
            linesWritten += 1
        podsWritten += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        print >> sys.stderr, linesRead, "lines read,", podsQueued, "pods queued,", podsDone, "pods filtered,", podsSorted, "pods sorted,", podsWritten, "pods written,", linesWritten, "good lines written."

def lineReader(fileObj):
    line = fileObj.readline()
    while len(line) >= 1:
        yield line
        line = fileObj.readline()

#########################################################################################################################


### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Input vcf file", action = "store", required = False)
parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--verbose", help="Verbose output.", action = "store_true")

parser.add_argument("-if", "--inputGenoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","alleles"])
parser.add_argument("-of", "--outputGenoFormat", help="Genotype format for output", action = "store", choices = ("phased","diplo","alleles","coded","count"), default = "phased")


#specific samples
parser.add_argument("-s", "--samples", help="sample names (separated by commas)", action='store')
#populations
parser.add_argument("-p", "--pop", help="Pop name and optionally sample names (separated by commas)", action='append', nargs="+", metavar=("popName","[samples]"))
parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
parser.add_argument("--keepAllSamples", help="Keep all samples (not just specified populations)", action='store_true')

#other
parser.add_argument("--haploid", help="Samples that are haploid (comma separated)", action = "store", metavar = "sample names")

#contigs
parser.add_argument("--include", help="include contigs (separated by commas)", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs (separated by commas)", action='store')
parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')


#further site filtering arguments
parser.add_argument("--minCalls", help="Minimum number of good genotype calls", type=int, action = "store", default = 1, metavar = "integer")
parser.add_argument("--minAlleles", help="Minimum number of alleles at a site", type=int, action = "store", default = 1, metavar = "integer")
parser.add_argument("--maxAlleles", help="Maximum number of alleles at a site", type=float, action = "store", default = "inf", metavar = "integer or 'inf'")
parser.add_argument("--minVarCount", help="Minimum number of instances for rare vaiants", type=int, action = "store", default = None, metavar = "integer")
parser.add_argument("--maxHet", help="Maximum proportion of heterozygous genotypes", type=float, action = "store", default = None, metavar = "proportion")
parser.add_argument("--minFreq", help="Minimum variant frequency", type=float, action = "store", default = None, metavar = "freqency")
parser.add_argument("--maxFreq", help="Maximum variant frequency", type=float, action = "store", default = None, metavar = "frequency")
parser.add_argument("--HWE", help="Hardy-Weinberg equalibrium test P-value and side", action = "store", nargs = 2, metavar = ("P-value", "'top'/'bottom'/'both'"))

#population-specific filtering arguments
parser.add_argument("--minPopCalls", help="Minimum number of good genotype calls per pop (comma separated)", action = "store", metavar = "integer")
parser.add_argument("--fixedDiffs", help="Only variants where differences are fixed between pops", action = "store_true")

#minimum distance for thinning
parser.add_argument("--thinDist", help="Allowed distance between sites for thinning", type=int, action = "store", metavar = "integer")

parser.add_argument("--skipChecks", help="Skip genotype checks to speed things up", action = "store_true")
parser.add_argument("--podSize", help="Lines to analyse in each thread simultaneously", type=int, action = "store", default = 100000)




args = parser.parse_args()

infile = args.infile
outfile = args.outfile

samples = args.samples
if samples: samples = samples.split(",")

include = []
exclude = []
if args.include: include += args.include.split(",")
if args.exclude: exclude += args.exclude.split(",")

if args.includeFile:
    with open(args.includeFile, 'r') as includeFile:
        include += includeFile.read().split()
        
if args.excludeFile:
    with open(args.excludeFile, 'r') as excludeFile:
        exclude += excludeFile.read().split()

if len(include) >= 1:
    include = set(include)
    print >> sys.stderr, "\nIncluding", len(include), "contigs.\n"
else: include = False

if len(exclude) >= 1:
    exclude = set(exclude)
    print >> sys.stderr, "\nExcluding", len(exclude), "contigs.\n"
else: exclude = False


minCalls = args.minCalls
minAlleles = args.minAlleles
maxAlleles = args.maxAlleles
minVarCount = args.minVarCount
maxHet = args.maxHet
minFreq = args.minFreq
maxFreq = args.maxFreq
fixed = args.fixedDiffs

if args.HWE:
    HWE_P = float(args.HWE[0])
    HWE_side = args.HWE[1]
else:
    HWE_P = HWE_side = None

popDict = {}
popNames = []
minPopCalls = None
if args.pop:
    
    for pop in args.pop:
        popNames.append(pop[0])
        popDict[pop[0]] = [] if len(pop)==1 else pop[1].split(",")
    
    if args.popsFile:
        with open(args.popsFile, "r") as pf: 
            for line in pf:
                ind,pop = line.split()
                if pop not in popDict: popDict[pop] = [ind]
                elif ind not in popDict[pop]: popDict[pop].append(ind)
    
    if args.minPopCalls: minPopCalls = dict(zip([pop[0] for pop in args.pop],
                                                [int(i) for i in args.minPopCalls.split(",")]))

nProcs = args.threads
verbose = args.verbose

##########################################################################################################################

### open files

if infile:
    if infile[-3:] == ".gz":
        In = gzip.open(infile, "r")
    else:
        In = open(infile, "r")
else:
    In = sys.stdin


if outfile:
    if outfile[-3:] == ".gz":
        Out = gzip.open(outfile, "w")
    else:
        Out = open(outfile, "w")
else:
    Out = sys.stdout


### read through header for all input files

headLine = In.readline()
headers = headLine.split()

#check specified samples are in first file. Otherwise use this entire set

allSamples = headers[2:]

if samples:
    for sample in samples:
        assert sample in allSamples, "Sample name not in header: " + sample
elif args.pop and not args.keepAllSamples:
    samples = sorted(list(set(popDict.values())))
else: samples = allSamples

for popName in popNames:
    for sample in popDict[popName]:
        assert sample in samples, "Specified population includes an unrecognised sample: " + sample

Out.write("\t".join(headers[0:2] + samples) + "\n")

ploidyDict = dict(zip(allSamples,[2]*len(allSamples)))

if args.haploid:
    for sample in args.haploid.split(","):
        ploidyDict[sample] = 1

##########################################################################################################################

#counting stat that will let keep track of how far we are
linesRead = 0
podsQueued = 0
podsDone = 0
podsSorted = 0
podsWritten = 0
linesWritten = 0


'''Create queues to hold the data. One will hold the pod info to be passed to the parser'''
inQueue = SimpleQueue()
#one will hold the results (in the order they come)
doneQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for parser. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s)
This one reads from the pod queue, passes each line some analysis function(s), gets the results and sends to the result queue'''
for x in range(nProcs):
    worker = Process(target=analysisWrapper,args=(inQueue,doneQueue,args.inputGenoFormat,args.outputGenoFormat,headers,include,exclude,samples,minCalls,minPopCalls,
                    minAlleles,maxAlleles,minVarCount,maxHet,minFreq,maxFreq,
                    HWE_P,HWE_side,popDict,ploidyDict,fixed,args.skipChecks,args.thinDist,))
    worker.daemon = True
    worker.start()

'''start two threads for sorting and writing the results'''
worker = Thread(target=sorter, args=(doneQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start one Process for sorting and writing the results'''
worker = Thread(target=writer, args=(writeQueue,Out,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()


########################################################################################################################

#place lines into pods
#pass pods on to processor(s)
podSize = args.podSize

lineGen = lineReader(In)
pod = []
podNumber = 0

for line in lineGen:
    linesRead += 1
    pod.append((linesRead, line))
    
    if linesRead % podSize == 0:
        inQueue.put((podNumber,pod))
        if verbose:
            print >> sys.stderr, "Pod", podNumber, "sent for analysis..."
        podNumber += 1
        podsQueued += 1
        pod = []


#run remaining lines in pod

if len(pod) > 0:
  inQueue.put((podNumber,pod))
  podsQueued += 1
  if verbose:
        print >> sys.stderr, "Pod", podNumber, "sent for analysis..."


#Wait for analysis to finish
while podsWritten < podsQueued:
    sleep(1)

sleep(10)

In.close()
if Out is not sys.stdout:
    Out.close()

