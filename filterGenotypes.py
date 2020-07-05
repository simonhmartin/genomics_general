#!/usr/bin/env python

import argparse, sys, gzip, string

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue

from threading import Thread

import genomics

from time import sleep


#######################################################################################################################

'''main worker function. This will watch the inQueue for pods, and pass lines from these pods to be parsed and filtered, before packaging back into a pod and sending on to the resultQueue'''
def analysisWrapper(inQueue,outQueue,inputGenoFormat,outputGenoFormat,alleleOrder,headers,include,exclude,samples,
                    minCalls,minPopCalls,minAlleles,maxAlleles,minPopAlleles,maxPopAlleles,minVarCount,maxHet,minFreq,maxFreq,
                    HWE_P,HWE_side,popDict,ploidyDict,fixed,nearlyFixedDiff,forcePloidy,partialToMissing,thinDist,noTest):
    sampleIndices = [headers.index(s) for s in samples]
    #dict of precomputed genotypes
    precompGTs = dict([(s, dict(),) for s in samples]) if precomp else None
    while True:
        podNumber,inPod = inQueue.get()
        if verbose: sys.stderr.write("Pod {} received for analysis.\n".format(podNumber))
        outPod = []
        lastScaf = None
        for lineData in inPod:
            lineNumber,line = lineData
            #if verbose: print >> sys.stderr, "Analysing line", lineNumber
            objects = line.split()
            if (include and objects[0] not in include) or (exclude and objects[0] in exclude): continue
            site = genomics.GenomeSite(genotypes=[objects[i] for i in sampleIndices], sampleNames=samples, popDict=popDict,
                                       ploidyDict=ploidyDict, genoFormat=inputGenoFormat, forcePloidy=forcePloidy, partialToMissing=partialToMissing, precompGTs=precompGTs)
            goodSite = True
            if thinDist:
                pos = int(objects[1])
                if lastScaf != objects[0]:
                    lastPos = pos
                    lastScaf = objects[0]
                    goodSite = False
                elif pos - lastPos < thinDist: goodSite = False
            if goodSite and not noTest: goodSite = genomics.siteTest(site,samples=samples,minCalls=minCalls,minPopCalls=minPopCalls,
                                               minAlleles=minAlleles,maxAlleles=maxAlleles,minPopAlleles=minPopAlleles,maxPopAlleles=maxPopAlleles,
                                               minVarCount=minVarCount,maxHet=maxHet,minFreq=minFreq,maxFreq=maxFreq,HWE_P=HWE_P,HWE_side=HWE_side,
                                               fixed=fixed,nearlyFixedDiff=nearlyFixedDiff)
            if goodSite:
                outLine = "\t".join(objects[:2] + [str(g) for g in site.asList(samples, mode=outputGenoFormat, alleleOrder=alleleOrder)]) + "\n"
                outPod.append((lineNumber,outLine))
                if thinDist: lastPos = int(objects[1])
            #if verbose: print >> sys.stderr, objects[0], objects[1], "passed: ", goodSite
        outQueue.put((podNumber,outPod))
        if verbose: sys.stderr.write("Pod {} analysed; sent to sorter.\n".format(podNumber))




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
            sys.stderr.write("Sorter received pod {}\n".format(podNumber))
        if podNumber == expect:
            writeQueue.put((podNumber,donePod))
            podsSorted += 1
            if verbose:
                sys.stderr.write("Pod {} sent to writer\n".format(podNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    donePod = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,donePod))
                    podsSorted += 1
                    if verbose:
                        sys.stderr.write("Pod {} sent to writer\n".format(expect))
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
            sys.stderr.write("Writer received pod {}.\n".format(podNumber))
        for thing in donePod:
            lineNumber,outLine = thing
            out.write(outLine)
            linesWritten += 1
        podsWritten += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} lines read | {} pods queued | {} pods filtered | {} pods sorted | {} pods written | {} good lines written.\n".format(linesRead, podsQueued, podsDone, podsSorted, podsWritten, linesWritten))


#########################################################################################################################


### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Input vcf file", action = "store", required = False)
parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--verbose", help="Verbose output.", action = "store_true")

parser.add_argument("-if", "--inputGenoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","alleles"], default="phased")

#ouput
parser.add_argument("-of", "--outputGenoFormat", action = "store", default = "phased",
                    choices = ("phased","diplo","bases","alleles","coded","count"), help="Genotype format for output")

parser.add_argument("--alleleOrder", action = "store", default = None,
                    choices = ("freq",), help="Order sample alleles by frequency when outputting 'bases' or 'alleles'")

#specific samples
parser.add_argument("-s", "--samples", help="sample names (separated by commas)", action='store')
parser.add_argument("--excludeSamples", help="sample names (separated by commas)", action='store')
#populations
parser.add_argument("-p", "--pop", help="Pop name and optionally sample names (separated by commas)", action='append', nargs="+", metavar=("popName","[samples]"))
parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
parser.add_argument("--keepAllSamples", help="Keep all samples (not just specified populations)", action='store_true')

#other
parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
parser.add_argument("--forcePloidy", help="Force genotypes to specified ploidy", action = "store_true")
parser.add_argument("--partialToMissing", help="Set partially missing genotypes to completely missing", action = "store_true")

#contigs
parser.add_argument("--include", help="include contigs", nargs = "+", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs", nargs = "+", action='store')
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
parser.add_argument("--minPopCalls", help="Minimum number of good genotype calls per pop", nargs="+", action = "store", type=int)
parser.add_argument("--minPopAlleles", help="Minimum number of alleles per site per pop", nargs="+", action = "store", type=int)
parser.add_argument("--maxPopAlleles", help="Maximum number of alleles per site per pop", nargs="+",action = "store", type=int)
parser.add_argument("--fixedDiffs", help="Only variants where differences are fixed between pops", action = "store_true")
parser.add_argument("--nearlyFixedDiff", help="Only variants where frequency diff between any pops is > x", action = "store", type=float)

#minimum distance for thinning
parser.add_argument("--thinDist", help="Allowed distance between sites for thinning", type=int, action = "store", metavar = "integer")

parser.add_argument("--podSize", help="Lines to analyse in each thread simultaneously", type=int, action = "store", default = 10000)
parser.add_argument("--noPrecomp", help="Do not use precomputed genotypes shortcut", action = "store_true")
parser.add_argument("--noTest", help="Output all lines (for debugging mostly)", action = "store_true")



args = parser.parse_args()

infile = args.infile
outfile = args.outfile

include = args.include if args.include else []
exclude = args.exclude if args.exclude else []

if args.includeFile:
    with open(args.includeFile, 'r') as includeFile:
        include += includeFile.read().split()
        
if args.excludeFile:
    with open(args.excludeFile, 'r') as excludeFile:
        exclude += excludeFile.read().split()

if len(include) >= 1:
    include = set(include)
    sys.stderr.write("\nIncluding {} contigs.\n".format(len(include)))
else: include = False

if len(exclude) >= 1:
    exclude = set(exclude)
    sys.stderr.write("\nExcluding {} contigs.\n".format(len(exclude)))
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
minPopCallsDict = None
minPopAllelesDict = None
maxPopAllelesDict = None
if args.pop:
    
    for pop in args.pop:
        popNames.append(pop[0])
        popDict[pop[0]] = [] if len(pop)==1 else pop[1].split(",")
    
    if args.popsFile:
        with open(args.popsFile, "rt") as pf: 
            for line in pf:
                ind,pop = line.split()
                if pop in popDict and ind not in popDict[pop]: popDict[pop].append(ind)
    
    if args.minPopCalls:
        minPopCalls = args.minPopCalls
        if len(minPopCalls) == 1: minPopCalls = minPopCalls*len(popNames)
        assert len(minPopCalls) == len(popNames)
        minPopCallsDict = dict(zip(popNames,minPopCalls))
    
    if args.minPopAlleles:
        minPopAlleles = args.minPopAlleles
        if len(minPopAlleles) == 1: minPopAlleles = minPopAlleles*len(popNames)
        assert len(minPopAlleles) == len(popNames)
        minPopAllelesDict = dict(zip(popNames,minPopAlleles))
        if args.maxPopAlleles == None: maxPopAllelesDict = dict(zip(popNames,[4]*len(popNames)))
    
    if args.maxPopAlleles:
        maxPopAlleles = args.maxPopAlleles
        if len(maxPopAlleles) == 1: maxPopAlleles = maxPopAlleles*len(popNames)
        assert len(maxPopAlleles) == len(popNames)
        maxPopAllelesDict = dict(zip(popNames,maxPopAlleles))
        if args.minPopAlleles == None: minPopAllelesDict = dict(zip(popNames,[0]*len(popNames)))
    
nProcs = args.threads
verbose = args.verbose

precomp = not args.noPrecomp

##########################################################################################################################

### open files

if infile:
    if infile[-3:] == ".gz":
        In = gzip.open(infile, "rt")
    else:
        In = open(infile, "rt")
else:
    In = sys.stdin


if outfile:
    if outfile[-3:] == ".gz":
        Out = gzip.open(outfile, "wt")
    else:
        Out = open(outfile, "wt")
else:
    Out = sys.stdout


### read through header for all input files

headLine = In.readline()
headers = headLine.split()

#check specified samples are in first file. Otherwise use this entire set

allSamples = headers[2:]

samples = args.samples.split(",") if args.samples else None
exSamples = args.excludeSamples.split(",") if args.excludeSamples else []

if samples is not None:
    for sample in samples:
        assert sample in allSamples, "Sample name not in header: " + sample
elif args.pop and not args.keepAllSamples:
    samples = [i for j in popDict.values() for i in j]
    assert len(set(samples)) == len(samples), "Populations cannot share the same sample"
else: samples = allSamples

samples = [s for s in samples if s not in exSamples]

if minCalls: assert minCalls <= len(samples), "Minimum calls is greater than number of specified samples."

for popName in popNames:
    popDict[popName] = [s for s in popDict[popName] if s not in exSamples]
    for sample in popDict[popName]:
        assert sample in allSamples, "Sample name not in header: " + sample

if args.ploidy is not None:
    ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(samples)
    assert len(ploidy) == len(samples), "Incorrect number of ploidy values supplied."
    ploidyDict = dict(zip(samples,ploidy))
elif args.ploidyFile is not None:
    with open(args.ploidyFile, "rt") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
else: ploidyDict = dict(zip(samples,[None]*len(samples)))

if args.outputGenoFormat != "bases":
    Out.write("\t".join(headers[0:2] + samples) + "\n")
else:
    assert args.ploidy != None, "Ploidy must be specified."
    outSamples = [sample + "_" + letter for sample in samples for letter in string.ascii_uppercase[:ploidyDict[sample]]]
    Out.write("\t".join(headers[0:2] + outSamples) + "\n")



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
    worker = Process(target=analysisWrapper,args=(inQueue,doneQueue,args.inputGenoFormat,args.outputGenoFormat,args.alleleOrder,
                                                  headers,include,exclude,samples,minCalls,minPopCallsDict,minAlleles,maxAlleles,
                                                  minPopAllelesDict,maxPopAllelesDict,minVarCount,maxHet,minFreq,maxFreq,
                                                  HWE_P,HWE_side,popDict,ploidyDict,fixed,args.nearlyFixedDiff,args.forcePloidy,
                                                  args.partialToMissing,args.thinDist,args.noTest,))
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

pod = []
podNumber = 0

for line in In:
    linesRead += 1
    pod.append((linesRead, line))
    
    if linesRead % podSize == 0:
        inQueue.put((podNumber,pod))
        if verbose:
            sys.stderr.write("Pod {} sent for analysis...\n".format(podNumber))
        podNumber += 1
        podsQueued += 1
        pod = []


#run remaining lines in pod

if len(pod) > 0:
  inQueue.put((podNumber,pod))
  podsQueued += 1
  if verbose:
        sys.stderr.write("Pod {} sent for analysis...\n".format(podNumber))


#Wait for analysis to finish
while podsWritten < podsQueued:
    sleep(1)

sleep(10)

In.close()
if Out is not sys.stdout:
    Out.close()

