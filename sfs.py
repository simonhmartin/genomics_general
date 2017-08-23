
'''Get frequency spectrum for each defined population, as well as the joint fs for all pairs and trips if requested'''

import sys,random,itertools,gzip,argparse
from time import sleep
import numpy as np
import genomics



ambigDict = {"K":["G","T"], "M":["A","C"], "R":["A","G"], "S":["C","G"], "W":["A","T"], "Y":["C","T"], "A":["A","A"], "C":["C","C"], "G":["G","G"], "T":["T","T"], "N":["N","N"] }

def haplo(calls):
    output = []
    for call in calls:
        if call in "ACGTNKRMSWY":
            output += ambigDict[call]
        else:
            print "WARNING", call, "is not recognised as a valid base or ambiguous base"
            output += ["N","N"]
    return output

def unique(things):
    output = list(set(things))
    output.sort()
    return output

def uniqueAlleles(bases):
    output = unique([i for i in bases if i in "ACGT"])
    return output


def initializeFS(nHapDict,popNames = None, doPairs = False, doTrios = False, nIntervals = 1):
    
    if not popNames:
        popNames = sorted(nHapDict.keys())
    
    ncat = [nHapDict[popName] + 1 for popName in popNames]
    ncat = dict(zip(popNames,ncat))
    
    FS = {}
    
    for popName in popNames:
        FS[popName] = np.zeros([ncat[popName]] + [nIntervals])

    if doPairs:
        pairs = list(itertools.combinations(popNames,2))
        for pair in pairs:
            FS["_".join(p for p in pair)] = np.zeros([ncat[p] for p in pair] + [nIntervals])

    if doTrios:
        trios = list(itertools.combinations(popNames,3))
        for trio in trios:
            FS["_".join(p for p in trio)] = np.zeros([ncat[p] for p in trio] + [nIntervals])
    
    return FS


def derCount(popBases, outBases):
    outAllele = uniqueAlleles(outBases)
    if len(outAllele) != 1:
        return None
    outAllele = outAllele[0]
    
    allPopBases = [b for pb in popBases for b in pb]
    popAlleles = uniqueAlleles(allPopBases)
    if len(popAlleles) > 2 or len(popAlleles) < 1:
        return None
    
    if len(popAlleles) == 1 and popAlleles[0] == outAllele[0]:
        return [0]*len(popBases)
    
    derAllele = [a for a in popAlleles if a != outAllele][0]
    
    derCounts = [p.count(derAllele) for p in popBases]
    
    return derCounts


def jitter(numbers, maxDist = 0.001):
  return [i + random.uniform(-maxDist, maxDist) for i in numbers]

def minorCount(popBases):
    allPopBases = [b for pb in popBases for b in pb]
    popAlleles = uniqueAlleles(allPopBases)
    
    if len(popAlleles) > 2:
        return None
    
    if len(popAlleles) == 1:
        return [0]*len(popBases)
    
    #we count up and add a tint jitter to ensure that one allele is always lower, even if there is a tie
    totalCounts = jitter([allPopBases.count(a) for a in popAlleles]) 
    
    minor = popAlleles[totalCounts.index(min(totalCounts))]
        
    minorCounts = [p.count(minor) for p in popBases]
    
    return minorCounts


def fsTableString(array):
    #first get number of dimensions
    dims = array.shape
    nPop = len(dims) - 1 # the last dimension is the number of intervals considered
    nIntervals = dims[-1]
    output = []
    #add header line
    #an allele freq column for each population, and number of sites for each interval
    output.append(",".join(["freq" for x in range(nPop)] + ["sites" for x in range(nIntervals)]))
    #now get the range for each population and use this to go through the array
    ranges = map(range,dims[:-1])
    #now we get all combos from itertools.
    #the * is a bit of python magic that feeds each range to itertools.product as a separate argument
    coords = itertools.product(*ranges)
    for coord in coords:
        #write the coords (repective freqs for each pop), and then the site counts for each interval
        output.append(",".join([str(x) for x in coord] + [str(int(array[tuple(list(coord) + [i])])) for i in range(nIntervals)]))
    
    return "\n".join(output)


def fsDadiString(array):

    #get dimensions
    dims = array.shape[:-1] # we ignore last dim, because there cannot be multiple intervals
    output = []
    #add dimensions
    output.append(" ".join([str(x) for x in dims]))
    #flatten out array into single list of length equal to the size of the array (this preserves the order of elements correctly for dadi (I THINK))
    siteCountList = list(array.reshape(array.size))
    #add to output
    output.append(" ".join([str(x) for x in siteCountList]))    
    return "\n".join(output)



def whichInterval(scaffold, position, scafIntervals, intervalPosDict):
    #this is very much a bespoke function for this script.
    #if your scaffold and position are in a defined interval or several intervals, this will return that interval number(s)
    try: scafInts = scafIntervals[scaffold]
    except: return None
    return [i for i in scafInts if position in intervalPosDict[i]]




'#####################################################################################################################'

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Data format for output", action = "store", choices = ("phased","diplo","alleles"), default = "phased")

parser.add_argument("-p", "--pop", help="Pop name and optionally sample names (separated by commas)", action='append', nargs="+", metavar=("popName","[samples]"))
parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")

parser.add_argument("--doPairs", help="Do pairwise fs for all pairs.", action="store_true")
parser.add_argument("--doTrios", help="Do trio fs for all trios.", action="store_true")
parser.add_argument("-S", "--subSample", help="Subsample number for each population", action='store', required = False, metavar=("number of individuals"))
parser.add_argument("--pref", help="Prefix for output files", action = "store", required = False, default = "")
parser.add_argument("--suff", help="Suffix for output files", action = "store", required = False, default = ".sfs")

parser.add_argument("--lastPopIsOutgroup", help="use last population as outgroup for polarizing unfloded spectrum", action='store_true')

parser.add_argument("--includeMono", help="Include counts at monomorphic sites (0 or all)", action="store_true")
parser.add_argument("--dadiFormat", help="Format output sfs for dadi", action="store_true")

#contigs
parser.add_argument("--include", help="include contigs (separated by commas)", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs (separated by commas)", action='store')
parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--intervalsFile", help="File of intervals, there will be a separate sfs for each", action='store')

parser.add_argument("-R", "--report", help="How often to report progress.", action='store', required = False, metavar=("number of lines"), default = 100000)
parser.add_argument("--verbose", help="Verbose output", action="store_true")
parser.add_argument("--test", help="Verbose output", action="store_true")


args = parser.parse_args()
#args = parser.parse_args("-c /home/simon/software/egglib_test/test.geno -p mel ros.CJ531,ros.CJ533,ros.CJ546,ros.CJ2071 -p cyd cyd.CJ553,cyd.CJ560,cyd.CJ564,cyd.CJ565 -O hec.JM273 --doPairs --pref test".split())

polarized = args.lastPopIsOutgroup
if not polarized: sys.stderr.write("\nNo outgroup provided. Minor allele frequency will be used.\n")

doPairs = args.doPairs

doTrios = args.doTrios

subSample = args.subSample
if subSample: ss=[int(x) for x in subSample.split(",")]

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

intervalsFile = args.intervalsFile

assert not (args.dadiFormat and intervalsFile), "You cannot specify separate intervale if outputting dadi format."

report = int(args.report)

verbose = args.verbose


'###########################################################################################################################'


### get files

if args.genoFile: genoFile = gzip.open(args.genoFile, "r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin

line = genoFile.readline()
header = line.split()

#IF INTERVALS ARE SPECIFIED, STORE THESE
if intervalsFile:
    print >> sys.stderr, "Making interval index..."
    with open(intervalsFile, "rU") as intFile:
        intervalsList = [line.split() for line in intFile.readlines()]
    nIntervals = len(intervalsList)
    #dicts giving scaffold and set of positions for all intervals
    intervalScafDict = dict([(i, intervalsList[i][0]) for i in xrange(nIntervals)])
    intervalPosDict = dict([(i, set(xrange(int(intervalsList[i][1]), int(intervalsList[i][2])+1))) for i in xrange(nIntervals)])
    #set of scaffolds so we can check whether to search further
    intervalScafs = set(intervalScafDict.values())
    #dict giving the intervals that a scaffold contains - speed up search
    scafIntervals = dict([(scaf, set([i for i in xrange(len(intervalScafDict)) if intervalScafDict[i] == scaf])) for scaf in intervalScafs])

else: nIntervals = 1

'###############################################################################################################'

#parse pops

popDict = {}
popNames = []
if args.pop:
    for pop in args.pop:
        popNames.append(pop[0])
        popDict[pop[0]] = [] if len(pop)==1 else pop[1].split(",")
    
    if args.popsFile:
        with open(args.popsFile, "r") as pf: 
            for line in pf:
                ind,pop = line.split()
                if pop in popDict and ind not in popDict[pop]: popDict[pop].append(ind)
    
    allSamples = [s for popName in popDict for s in popDict[popName]]
else:
    popNames = ["all"]
    popDict = {"all": header[2:]}


for popName in popNames: assert len(popDict[popName]) >= 1, "Population {} has no samples".format(popName) 

allSamples = [s for popName in popDict for s in popDict[popName]]

if polarized:
    inPopNames = popNames[:-1]
    inPopDict = popDict.copy()
    outgroupName = popNames.pop(-1)
    OGsamples = popDict.pop(outgroupName) #haha don't get confused. Pop is also a function:)
else:
    inPopNames = popNames
    inPopDict = popDict

sampleIndices = [header.index(s) for s in allSamples]

if args.ploidy is not None:
    ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allSamples)
    assert len(ploidy) == len(allSamples), "Incorrect number of ploidy values supplied."
    ploidyDict = dict(zip(allSamples,ploidy))
elif args.ploidyFile is not None:
    with open(args.ploidyFile, "r") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
else: ploidyDict = dict(zip(allSamples,[2]*len(allSamples)))


if args.verbose:
    for popName in inPopNames: sys.stderr.write("\n"+ popName + ":\n" + ",".join(inPopDict[popName]) + "\n")
    if polarized: sys.stderr.write("\noutgroup = "+ outgroupName + ":\n" + ",".join(OGsamples) + "\n")


P = len(inPopNames)

if doPairs:
    pairs = list(itertools.combinations(inPopNames,2))
    pairNames = ["_".join(p for p in pair) for pair in pairs]
    pairs = dict(zip(pairNames,pairs))

if doTrios:
    trios = list(itertools.combinations(inPopNames,3))
    trioNames = ["_".join(p for p in trio) for trio in trios]
    trios = dict(zip(trioNames,trios))


#samples per pop

nHapDict = dict([(popName, sum([ploidyDict[sample] for sample in popDict[popName]])) for name in inPopNames]) 

#if sub sampling, first check that all pops have at least that number of samples
if subSample is not None:
    if len(subSample) == 1: ss = ss*len(inPopNames)
    else: assert len(ss) == len(inPopNames), "subsample list must match number of ingroup populations"
    ssDict = dict(zip(inPopNames, ss))
    for p in inPopNames:
        assert nHapDict[p] >= ssDict[p], "Population {} has fewer than {} haplotypes ({}).".format(p,ssDict[p],str(nHapDict[p])) 
    nHapDict = ssDict

#initialize arrays of allele counts
arrays = initializeFS(nHapDict, inPopNames, doPairs, doTrios, nIntervals)


linesDone = 0
sitesAnalysed = 0
snpsCounted = 0

### for each line, check if its a scaf we want
for line in genoFile:
    
    linesDone += 1
    if linesDone % report == 0: print >> sys.stderr, linesDone, "lines done..."
    
    objects = line.split()
    
    if (include and objects[0] not in include) or (exclude and objects[0] in exclude): continue
    
    #if there are intervals, check whether the site matches any
    
    if intervalsFile: siteIntervals = whichInterval(objects[0], int(objects[1]), scafIntervals, intervalPosDict)
    else: siteIntervals = [0]
    
    if siteIntervals:
                
        site = genomics.GenomeSite(genotypes=[objects[i] for i in sampleIndices], sampleNames=allSamples, popDict=popDict,
                                ploidyDict=ploidyDict, genoFormat=args.genoFormat)
            
        popBases = dict([(popName, [b for b in site.asList(samples=popDict[popName], mode="bases") if b in "ACGT"],) for popName in inPopNames])
        
        if args.test:
            sleep(0.25)
            for popName in inPopNames:
                sys.stderr.write("\n\n" + popName + ": " + " ".join(popBases[popName]))
        
        goodData = True
        
        #subsample if necessary
        if subSample:
            for popName in popNames:
                try: popBases[popName] = random.sample(popBases[popName],ssDict[popName])
                except: goodData = False
            #if we get an error its because there's insufficient data in at least one pop, so we ignore the whole site
            # in theory we could modify this part to use info for the pops it can - might be necessary when sites are limited
        #if we didn't subsample, check that we have the right number of calls from each pop
        else:
            try: assert np.all([len(popBases[popName]) == nHapDict[popName] for popName in inPopNames])
            except: goodData = False
            # as in the subsample above, here we conservatively reject the entire line if one pop is missing data
            # we might want to change this in future
        
        if not goodData: continue
        
        #allele counts per pop
        countsArray = np.array([genomics.baseFreqs(popBases[popName], asCounts=True) for popName in inPopNames])
        
        #total allele counts summed over pops
        inCounts = np.sum(countsArray, axis = 0)
                
        #check mono or biallelic
        alleles = inCounts != 0
        if not 1 <= np.sum(alleles) <= 2: continue
        
        #if we get here, the site has sufficient data to be considered for analysis
        sitesAnalysed += 1
        
        # if there's an outgroup to polariize - we get the outgroup bases and then the derived counts, otherwise get minor allele counts
        #the allele we count (minor or derived) will be called the target
        if polarized:
            outBases = [b for b in site.asList(samples=OGsamples, mode="bases") if b in "ACGT"]
            
            if args.test: sys.stderr.write("\n"+outgroupName + ": " + " ".join(outBases))

            outCounts = genomics.baseFreqs(outBases, asCounts=True)
            outAlleles = outCounts != 0
            #check there is one outgroup allele and it is one of the ingroup alleles
            if sum(outAlleles) != 1 or sum(outAlleles & alleles) > 2: continue
            try: target = np.where(~outAlleles & alleles)[0][0]
            except: target = np.where(~alleles)[0][0] #error here means invariant site, so take any absent allele as target 
        
        else: target = np.where(inCounts == sorted(inCounts)[-2])[0][0]
        
        if args.test: sys.stderr.write("\nTarget base: {} ".format(target))
        
        popTargetCounts = dict(zip(inPopNames, countsArray[:,target]))
        
        if args.test:
            for popName in inPopNames:
                sys.stderr.write(popName + ": " + str(popTargetCounts[popName]) + " ")
        
        
        #if we get here we are going to add this data tothe SFS
        snpsCounted += 1
        
        for si in siteIntervals:
            for popName in inPopNames:
                arrays[popName][popTargetCounts[popName],si] += 1
            
            #for pairs
            if doPairs:
                for pairName in pairNames:
                    arrays[pairName][popTargetCounts[pairs[pairName][0]],popTargetCounts[pairs[pairName][1]], si] += 1
            
            #check for trios
            if doTrios:
                for trioName in trioNames:
                    arrays[trioName][popTargetCounts[trios[trioName][0]],popTargetCounts[trios[trioName][1]],popTargetCounts[trios[trioName][2]], si] += 1


if not args.includeMono:
    #if we're not including monomorphic sites, we must remove the zero class and the ns class where all samples have a derived allele
        
    for popName in inPopNames:
        arrays[popName][0,:] = arrays[popName][nHapDict[popName]] = np.NAN

    '''We do the same for pair spectra. Here only deleting the [0,0] and [ns,ns] classes.
        And then also the same for trios.'''

    if doPairs:
        for pairName in pairNames:
            arrays[pairName][0,0,:] = arrays[pairName][nHapDict[pairs[pairName][0]],nHapDict[pairs[pairName][1]]] = np.NAN

    if doTrios:
        for trioName in trioNames:
            arrays[trioName][0,0,0,:] = arrays[trioName][nHapDict[trios[trioName][0]],nHapDict[trios[trioName][1]],nHapDict[trios[trioName][2]]] = np.NAN


#write output files. A separate file for each fs
for name in arrays.keys():
    with open(args.pref + name + args.suff, "w") as out:
        if args.dadiFormat:
            out.write("#Sites analysed: " + str(sitesAnalysed) + "\n" + "#SNPs counted: " + str(snpsCounted) + "\n")
            out.write(fsDadiString(arrays[name]) + "\n")
        else: out.write(fsTableString(arrays[name]) + "\n")
