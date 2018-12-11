#!/usr/bin/env python

import argparse, sys, gzip, string

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread

import numpy as np

from scipy.stats import fisher_exact,chisquare

import genomics

from time import sleep


#######################################################################################################################

'''main worker function. This will watch the inQueue for pods, and pass lines from these pods to be parsed and filtered, before packaging back into a pod and sending on to the resultQueue'''
def analysisWrapper(inQueue,outQueue,inputGenoFormat,headers,include,exclude,group1inds,group2inds,permutations, permutationMaxP):
    
    samples = group1inds + group2inds
    sampleIndices = [headers.index(s) for s in samples]
    
    group1 = np.array([True]*len(group1inds) + [False]*len(group2inds))
    group2 = ~group1
    
    while True:
        podNumber,inPod = inQueue.get()
        if verbose: sys.stderr.write("Pod {} received for analysis...\n".format(podNumber))

        
        outPod = []
        for lineData in inPod:
            lineNumber,line = lineData
            #if verbose: print >> sys.stderr, "Analysing line", lineNumber
            objects = line.split()
            if (include and objects[0] not in include) or (exclude and objects[0] in exclude): continue
            site = genomics.GenomeSite(genotypes=[objects[i] for i in sampleIndices],
                                       sampleNames=samples,genoFormat=inputGenoFormat)
            
            alleles = site.alleles()
            
            if len(alleles) == 2:
                
                minorCount = np.array(site.asList(mode="count", countAllele = alleles[1], missing=-1))
                majorCount = np.array(site.asList(mode="count", countAllele = alleles[0], missing=-1))
                
                #get index for good genotypes and filter all by that
                idx = np.where(minorCount >= 0)[0]
                
                _group1_ = group1[idx]
                _group2_ = group2[idx]
                
                minorPresent = minorCount[idx] >= 1
                minorAbsent = ~minorPresent
                
                majorPresent = majorCount[idx] >= 1
                majorAbsent = ~majorPresent
                
                minorTable = np.array([[(minorPresent & _group1_).sum(),(minorAbsent & _group1_).sum()],
                                         [(minorPresent & _group2_).sum(),(minorAbsent & _group2_).sum()]])
                
                majorTable = np.array([[(majorPresent & _group1_).sum(),(majorAbsent & _group1_).sum()],
                                        [(majorPresent & _group2_).sum(),(majorAbsent & _group2_).sum()]])
                
                p_values = (fisher_exact(minorTable)[1], fisher_exact(majorTable)[1],)
                
                result = [min(p_values)]
                
                if permutations >= 1:
                    if permutationMaxP is None or result[0] <= permutationMaxP:
                        
                        table = minorTable if p_values[0] <= p_values[1] else majorTable 
                        
                        phi = chisquare(table, axis=None)[0]/table.sum()
                        
                        phi_permuted = []
                        for i in range(permutations):
                            newGroup1 = np.random.permutation(_group1_)
                            newGroup2 = ~newGroup1
                            
                            newTable = np.array([[(minorPresent & newGroup1).sum(),(minorAbsent & newGroup1).sum()],
                                                 [(minorPresent & newGroup2).sum(),(minorAbsent & newGroup2).sum()]])
                            
                            phi_permuted.append(chisquare(newTable, axis=None)[0]/table.sum())
                        
                        p_emp = (len([_phi_ for _phi_ in phi_permuted if _phi_ >= phi]) + 1.) / (permutations + 1.) 
                    
                    else: p_emp = np.NaN
                    
                    result.append(p_emp)
                
            elif permutations >= 1: result = [np.NaN]*2
            
            else: result = [np.NaN]
            
            outLine = "\t".join(objects[:2] + [str(round(x, 5)) for x in result]) + "\n"
            
            outPod.append((lineNumber,outLine))
        
        outQueue.put((podNumber,outPod))
        if verbose: sys.stderr.write("Pod {} analysed, sent to sorter.\n".format(podNumber))



'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, verbose):
    global podsDone
    global podsSorted
    sortBuffer = {}
    expect = 0
    while True:
        podNumber, donePod = doneQueue.get()
        podsDone += 1
        if verbose: sys.stderr.write("Sorter received pod {}\n".format(podNumber))

        if podNumber == expect:
            writeQueue.put((podNumber,donePod))
            podsSorted += 1
            if verbose: sys.stderr.write("Pod {} sent to writer...\n".format(podNumber))

            expect +=1
            #now check buffer for further results
            while True:
                try:
                    donePod = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,donePod))
                    podsSorted += 1
                    if verbose: sys.stderr.write("Pod {} sent to writer...\n".format(podNumber))
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
        if verbose: sys.stderr.write("Writer received pod {}\n".format(podNumber))
        for thing in donePod:
            lineNumber,outLine = thing
            out.write(outLine)
            linesWritten += 1
        podsWritten += 1


'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} lines read, {} pods queued, {} pods analysed, "\
                         "{} pods sorted, {} pods written, {} lines written.\n".format(linesRead, podsQueued,podsDone,
                                                                                      podsSorted,podsWritten,linesWritten))

#########################################################################################################################


### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Input vcf file", action = "store", required = False)
parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("--verbose", help="Verbose output.", action = "store_true")

parser.add_argument("-if", "--inputGenoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","alleles"], default="phased")

#phenotypes
parser.add_argument("--phenoFile", help="File of sample names and phenotypes", action = "store", required = True)
parser.add_argument("--phenoColumn", help="Column to use in phenoFile", action = "store", default=2)
parser.add_argument("--phenotypes", help="phenotype names to compare, if there are more than 2", action='store', nargs=2)

#analysis
parser.add_argument("--permutations", help="Permutations for empirical p value", action = "store", type=int, default=0)

parser.add_argument("--permutationMaxP", action = "store", type=float,
                    help="Only perform permutation if p is less or equal to this")

#contigs
parser.add_argument("--include", help="include contigs", nargs = "+", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs", nargs = "+", action='store')
parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')

parser.add_argument("--podSize", help="Lines to analyse in each thread simultaneously", type=int, action = "store", default = 100000)



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
    print >> sys.stderr, "\nIncluding", len(include), "contigs.\n"
else: include = False

if len(exclude) >= 1:
    exclude = set(exclude)
    print >> sys.stderr, "\nExcluding", len(exclude), "contigs.\n"
else: exclude = False


#get ind names and phenotypes

indNames = []
indPhenos = []

with open(args.phenoFile, "r") as pf:
    if args.phenoColumn:
        try: column = int(args.phenoColumn)
        except:
            pf_headers = pf.readline().split()
            assert args.phenoColumn in pf_headers, "\nPhenotype column name not recognised.\n"
            column = pf_headers.index(args.phenoColumn)
    else: column=1
    
    for line in pf:
        elements = line.split()
        indNames.append(elements[0])
        indPhenos.append(elements[column])

if not args.phenotypes:
    phenotypes=np.unique(indPhenos)
    assert len(phenotypes) == 2, "\nPlease specify two phenotypes to consider using --phenotypes\n"

else: assert args.phenotypes[0] in indPhenos and args.phenotypes[1] in indPhenos, "\nGiven phenotype not in phenotypes file.\t"

group1inds = [indNames[x] for x in range(len(indNames)) if indPhenos[x] == phenotypes[0]]
group2inds = [indNames[x] for x in range(len(indNames)) if indPhenos[x] == phenotypes[1]]


sys.stderr.write("\nGroup 1 ({}):\n {}\n".format(phenotypes[0], " ".join(group1inds)))
sys.stderr.write("\nGroup 2 ({}):\n {}\n".format(phenotypes[1], " ".join(group2inds)))

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

#check specified individuals are in first file. Otherwise use this entire set

allInds = headers[2:]

for ind in group1inds + group2inds:
    assert ind in allInds, "\nIndividual {} not found in file header\n".format(ind)

### write output header

if args.permutations >= 1:
    resultHeaders = ["p_fisher", "p_emp"]
    sys.stderr.write("\nEmpirical p values to be calculated using {} permutations\n\n".format(args.permutations))

else: resultHeaders = ["p_fisher"]

Out.write("\t".join(headers[:2] + resultHeaders) + "\n")

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
    worker = Process(target=analysisWrapper,args=(inQueue,doneQueue,args.inputGenoFormat,
                                                  headers,include,exclude,group1inds,group2inds,args.permutations, args.permutationMaxP))
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
        #pause here if we're waiting for too many pods to be analysed, sorted and written
        while podsDone - podsWritten >= args.threads*2: sleep(5)
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

