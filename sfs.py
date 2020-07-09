#!/usr/bin/env python

'''Get frequency spectrum for each defined population, as well as the joint fs for pairs and trips if requested'''

import sys,random,itertools,gzip,argparse
from time import sleep
import numpy as np
from collections import defaultdict
import genomics
import time


#function to take genotype data for a site and get the four base counts for each individual, separated by population
def getPopIndBaseCounts(siteData, genoFormat, allSamples, popDict, ploidyDict):
    site = genomics.GenomeSite(genotypes=[siteData["GTs"][name] for name in allSamples], sampleNames=allSamples,
                               popDict=popDict, ploidyDict=ploidyDict, genoFormat=args.genoFormat)
    
    popIndBaseCounts = dict([(popName, np.array([site.genotypes[indName].asBaseCounts() for indName in popDict[popName]]),) for popName in popNames])
    
    return popIndBaseCounts
    

def downSampleBaseCounts(baseCounts, N):
    return np.bincount(np.random.choice(np.repeat(np.arange(4),baseCounts),N,replace=False),minlength=4)

### NOTE ### below is some code to do the downsampling using the hypergeometric distribution
### however, it SOMETIMES fails, giving a different number of alleles to the number requested
### It is not detectably slower than the version above for moderate sample sizes
### I'm only leaving this here as a reminder to myself NOT to use it

#def downSampleBaseCounts(baseCounts, N):
    #return np.random.hypergeometric(baseCounts,baseCounts.sum()-baseCounts, N)


#function 
def getPopBaseCounts(popIndBaseCounts, popNames, subsampleDict=None, subsampleIndividuals=False):
    
    # get population base counts, and do the subsampling if necessary
    # This is currently conservative. If any one of the populations lacks sufficient good genotypes it will break
    # in theory we could modify this part to use info for the pops it can - might be necessary when sites are limited

    if subsampleDict:
        #this is subsampling by individual, and it is usually not necessary, but could beof interest
        if subsampleIndividuals:
            #index individuals with non-zero base counts
            popGoodInds = dict([(popName, np.where(popIndBaseCounts[popName].sum(axis=1) != 0)[0],) for popName in popNames])
            
            try: popBaseCountsArray = np.array([popIndBaseCounts[popName][random.sample(popGoodInds[popName],subsampleDict[popName]),:].sum(axis = 0) for popName in popNames])
            except: return
        else:
            try: popBaseCountsArray = np.array([downSampleBaseCounts(popIndBaseCounts[popName].sum(axis = 0), subsampleDict[popName]) for popName in popNames])
            except: return
    else:
        popBaseCountsArray = np.array([popIndBaseCounts[popName].sum(axis = 0) for popName in popNames])
    
    return popBaseCountsArray

#function to get the count of a target allele given an aray of base counts for each pop
#if outgroup base counts are given, then the derived allele is used, if possible
def getTargetCounts(popBaseCountsArray, outgroupBaseCounts=None, outgroupMono=True):
    #total allele counts summed over pops
    totalBaseCounts = popBaseCountsArray.sum(axis = 0)
    
    alleles = totalBaseCounts > 0
    
    if outgroupBaseCounts is not None:
        outAlleles = outgroupBaseCounts > 0
        allAlleles = alleles | outAlleles
    else: allAlleles = alleles
    
    #check mono or biallelic
    if not 1 <= allAlleles.sum() <= 2: return
    
    if outgroupBaseCounts is not None:
        nOutAlleles = outAlleles.sum()
        #check there is one outgroup allele and it is one of the ingroup alleles
        if nOutAlleles == 0 or (outgroupMono & nOutAlleles != 1): return
        #set target
        try: target = np.where(~outAlleles & alleles)[0][0]
        except: target = np.where(~alleles)[0][0] #error here means invariant site, so take first absent allele as target 
    
    #otherwise take the minor (second most frequent) allele as target
    else: target = totalBaseCounts.argsort()[-2]
    
    return popBaseCountsArray[:,target]


### Below you will see perhaps the cleverest piece of Python I will ever write.
### It creates either a single defaultdict that defaults to zeros, for the case of a 1D frequency spectrum
### Or nested dicts according to the number of dimensions specified.
### These are sparse spectra, because they only contain the values given to them.
### using defaultdicts means that the values (and their keys at each dimension) are automatically generated when being set
### Even though it works, I must admit that I onlly understand about 80% why ;)
class SparseFS(defaultdict):
    def __init__(self, dimensions=1, intervals=1):
        self.dimensions = dimensions
        self.intervals = intervals
        if dimensions == 1: super().__init__(lambda: np.zeros(intervals, dtype=int))
        else:
            super().__init__(lambda: SparseFS(dimensions-1, intervals))
    
    def getCount(self, freqs):
        if self.dimensions == 1: return self[freqs[0]]
        else: return self[freqs[0]].getCount(freqs[1:])
    
    def setCount(self, freqs, value=1):
        if self.dimensions == 1:
            self[freqs[0]] = value
        else: return self[freqs[0]].setCount(freqs[1:], value)
    
    def add(self, freqs, value=1):
        if self.dimensions == 1:
            self[freqs[0]] += value
        else: return self[freqs[0]].add(freqs[1:], value)
    
    def asChains(self, chains = list(), chain=list()):
        if self.dimensions == 1:
            for key in self.keys():
                yield chain + [key] + list(self[key])
        else:
            for key in self.keys():
               yield from self[key].asChains(chain=chain + [key])

#def fsTableString(array):
    ##first get number of dimensions
    #dims = array.shape
    #nPop = len(dims) - 1 # the last dimension is the number of intervals considered
    #nIntervals = dims[-1]
    #output = []
    ##add header line
    ##an allele freq column for each population, and number of sites for each interval
    #output.append("\t".join(["freq" for x in range(nPop)] + ["sites" for x in range(nIntervals)]))
    ##now get the range for each population and use this to go through the array
    #ranges = map(range,dims[:-1])
    ##now we get all combos from itertools.
    ##the * is a bit of python magic that feeds each range to itertools.product as a separate argument
    #coords = itertools.product(*ranges)
    #for coord in coords:
        ##write the coords (repective freqs for each pop), and then the site counts for each interval
        #output.append(",".join([str(x) for x in coord] + [str(int(array[tuple(list(coord) + [i])])) for i in range(nIntervals)]))
    
    #return "\n".join(output)


#def fsDadiString(array):

    ##get dimensions
    #dims = array.shape[:-1] # we ignore last dim, because there cannot be multiple intervals
    #output = []
    ##add dimensions
    #output.append(" ".join([str(x) for x in dims]))
    ##flatten out array into single list of length equal to the size of the array (this preserves the order of elements correctly for dadi (I THINK))
    #siteCountList = list(array.reshape(array.size))
    ##add to output
    #output.append(" ".join([str(x) for x in siteCountList]))    
    #return "\n".join(output)


'#####################################################################################################################'

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", help="Input file", action = "store")

parser.add_argument("--inputType", action = "store", choices = ("genotypes","baseCounts","targetCounts"),
                    default = "targetCounts", help="Data type of input")

parser.add_argument("--scafCol", action='store', type=int, default = 0,
                    help="Specify which column has scaffold information. -1 if no scaffold column")

parser.add_argument("--posCol", action='store', type=int, default = 1,
                    help="Specify which column has position information. -1 if no position column")

parser.add_argument("--firstSampleCol", action='store', type=int, default = 2,
                    help="Specify which column has information for the first sample/population")

parser.add_argument("--header", help="Input file header if none present (must match sample/population names)", action="store")

parser.add_argument("--genoFormat", action = "store", choices = ("phased","diplo","alleles"), default = "phased",
                    help="If input is genotypes - what format are they in")

parser.add_argument("-p", "--pop", action='append', nargs="+", metavar=("popName","[samples]"),
                    help="Pop name and optionally sample names (separated by commas)")
parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")

parser.add_argument("--FSpops", help="Pop or pops to make a sfs for", action = "append", type=str, nargs="+")
parser.add_argument("--doPairs", help="Do pairwise fs for all pairs.", action="store_true")
parser.add_argument("--doTrios", help="Do trio fs for all trios.", action="store_true")
parser.add_argument("--doQuartets", help="Do 4D fs for all quartets.", action="store_true")

parser.add_argument("--subsample", action='store', required = False,  nargs = "+", type=int,
                    metavar=("number of bases"), help="#haplotypes (or individuals) to subsample from each population")
parser.add_argument("--subsampleIndividuals", action='store_true',
                    help="Subsample by individual rather than treating haplotypes separately")

parser.add_argument("--pref", help="Prefix for output files", action = "store", required = False, default = "")
parser.add_argument("--suff", help="Suffix for output files", action = "store", required = False, default = ".sfs")
parser.add_argument("--pipe", help="Write all outputs to stdout", action = "store_true")

parser.add_argument("--polarized", help="use last population as outgroup for polarizing unfloded spectrum", action='store_true')
parser.add_argument("--outgroup", help="Outgroup for polarizing frequencies", action = "store", required = False)

#parser.add_argument("--dadiFormat", help="Format output sfs for dadi", action="store_true")

#contigs
parser.add_argument("--include", help="include contigs", nargs = "+", action='store')
parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
parser.add_argument("--exclude", help="exclude contigs", nargs = "+", action='store')
parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')

parser.add_argument("--regions", help="Regions defined as chrom or chrom:start-end. Will create a separate sfs for each", nargs = "+", action='store')
parser.add_argument("--regionsFile", help="File of regions with chrom and optionally start and end", action='store')

parser.add_argument("-R", "--report", help="How often to report progress.", action='store', required = False, metavar=("number of lines"), default = 100000)
parser.add_argument("--verbose", help="Verbose output", action="store_true")
parser.add_argument("--seed", action='store', type=int, help="Set random seed for subsampling", default=42)


args = parser.parse_args()

if not (args.polarized | args.outgroup) and args.inputType != "targetCounts":
    sys.stderr.write("\nNo outgroup provided. Minor allele frequency will be used.\n")

subsample = args.subsample

include = args.include if args.include else []
exclude = args.exclude if args.exclude else []

if args.includeFile:
    with open(args.includeFile, 'rt') as includeFile:
        include += includeFile.read().split()
        
if args.excludeFile:
    with open(args.excludeFile, 'rt') as excludeFile:
        exclude += excludeFile.read().split()

if len(include) >= 1:
    include = set(include)
    sys.stderr.write("\nIncluding {} contigs.\n".format(len(include)))
else: include = False

if len(exclude) >= 1:
    exclude = set(exclude)
    sys.stderr.write("\nExcluding {} contigs.\n".format(len(exclude)))
else: exclude = False

#assert not (args.dadiFormat and (args.regionsFile or args.regions)), "You cannot specify separate intervals if outputting dadi format."

report = int(args.report)

verbose = args.verbose

np.random.seed(args.seed)

'###########################################################################################################################'


### get files

if args.inputFile: inputFile = gzip.open(args.inputFile, "rt") if args.inputFile.endswith(".gz") else open(args.inputFile, "rt")
else: inputFile = sys.stdin


#IF INTERVALS ARE SPECIFIED, STORE THESE
if args.regions or args.regionsFile:
    if args.regions:
        intervals = genomics.Intervals(regions=args.regions)
    else:
        with open(args.regionsFile, "rt") as intFile:
            intervals = genomics.Intervals(tuples=[line.split() for line in intFile])
    nIntervals = len(intervals.chroms)
    sys.stderr.write("Recording SFS for {} intervals\n".format(nIntervals))
else:
    intervals = None
    nIntervals = 1

'###############################################################################################################'

#parse pop and individual data, if necessary
if args.inputType == "genotypes":
    genoFileReader = genomics.GenoFileReader(inputFile, headerLine=args.header,
                                             scafCol=args.scafCol, posCol=args.posCol, firstSampleCol=args.firstSampleCol)
    popDict = {}
    popNames = []
    if args.pop or args.FSpops:
        if args.pop:
            for pop in args.pop:
                popNames.append(pop[0])
                popDict[pop[0]] = [] if len(pop)==1 else pop[1].split(",")
        
        if args.FSpops:
            for pop in [p for pops in args.FSpops for p in pops]:
                if pop not in popNames:
                    popNames.append(pop)
                    popDict[pop] = []
        
        if args.popsFile:
            with open(args.popsFile, "r") as pf: 
                for line in pf:
                    ind,pop = line.split()
                    if pop in popDict and ind not in popDict[pop]: popDict[pop].append(ind)
        
        allSamples = [s for popName in popDict for s in popDict[popName]]
    else:
        #populations not specified, assume all individuals in one pop
        popNames = ["all"]
        popDict = {"all": genoFileReader.names}
    
    for popName in popNames: assert len(popDict[popName]) >= 1, "Population {} has no samples".format(popName) 
    
    allSamples = [s for popName in popDict for s in popDict[popName]]
    
    if args.ploidy is not None:
        ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allSamples)
        assert len(ploidy) == len(allSamples), "Incorrect number of ploidy values supplied."
        ploidyDict = dict(zip(allSamples,ploidy))
    
    elif args.ploidyFile is not None:
        with open(args.ploidyFile, "r") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
    
    else: ploidyDict = dict(zip(allSamples,[2]*len(allSamples)))
    
    #number of haplotypes per pop, determines size of the sfs
    nHapDict = dict([(popName, sum([ploidyDict[sample] for sample in popDict[popName]])) for popName in popNames]) 
    
    if args.verbose:
        for popName in popNames: sys.stderr.write("\n"+ popName + ":\n" + ",".join(popDict[popName]) + "\n")

elif args.inputType == "baseCounts":
    #assume input file contains comma-separated counts of each base
    genoFileReader = genomics.GenoFileReader(inputFile, headerLine=args.header,
                                             scafCol=args.scafCol, posCol=args.posCol, firstSampleCol=args.firstSampleCol)
    if args.pop or args.FSpops:
        
        popNames = []
        
        if args.pop:
            for pop in args.pop: popNames.append(pop[0])
        
        if args.FSpops:
            for pop in [p for pops in args.FSpops for p in pops]:
                if pop not in popNames: popNames.append(pop)
    else:
        popNames = genoFileReader.names

else:
    #otherwise we assume only the target base is given in the input file
    genoFileReader = genomics.GenoFileReader(inputFile, headerLine=args.header,
                                             scafCol=args.scafCol, posCol=args.posCol, firstSampleCol=args.firstSampleCol, type=int)
    if args.pop or args.FSpops:
        
        popNames = []
        
        if args.pop:
            for pop in args.pop: popNames.append(pop[0])
        
        if args.FSpops:
            for pop in [p for pops in args.FSpops for p in pops]:
                if pop not in popNames: popNames.append(pop)
    else:
        popNames = genoFileReader.names


sys.stderr.write("\nPopulations:\n")
sys.stderr.write(" ".join(popNames) + "\n")

#if polarizing, assume last population is outgroup
if (args.inputType == "genotypes" or args.inputType == "baseCounts") and (args.polarized or args.outgroup):
    outgroup = args.outgroup if args.outgroup else popNames[-1]
    inPopNames = [popName for popName in popNames if popName != outgroup]
    sys.stderr.write("\nFrequencies will be polarized assuming outgroup is {}\n".format(outgroup))
else:
    #otherwise just consider all is ingroups
    inPopNames = popNames


#if sub sampling
if subsample is not None:
    
    if len(subsample) == 1: subsample = subsample*len(inPopNames)
    else: assert len(subsample) == len(inPopNames), "subsample list ({}) must match number of ingroup populations ({}).".format(len(subsample),len(inPopNames))
    
    subsampleDict = dict(zip(inPopNames, subsample))
    
    if args.inputType == "genotypes":
        if not args.subsampleIndividuals:
            #assume that subsampleDict reflects number of haplotypes to subsample
            #check that all pops have at least that number of samples
            for p in inPopNames:
                assert nHapDict[p] >= subsampleDict[p], "Population {} has fewer than {} haplotypes ({}).".format(p,subsampleDict[p],str(nHapDict[p])) 
            nHapDict = subsampleDict
        else:
            globalPloidy = tuple(set([ploidyDict[ind] for pop in inPopNames for ind in popDict[pop]]))
            assert(len(globalPloidy) == 1), "Subsampling by individuals not possible with variable ploidy"
            #assume that subsampleDict reflects number of INDIVIDUALS to subsample
            _nHapDict_ = dict(zip(inPopNames, [s*globalPloidy[0] for s in subsample]))
            for p in inPopNames:
                assert nHapDict[p] >= _nHapDict_[p], "Population {} has fewer than {} haplotypes ({}).".format(p,_nHapDict_[p],str(nHapDict[p])) 
            nHapDict = _nHapDict_

else: subsampleDict = None

if args.inputType == "genotypes":
    nHapArray = np.array([nHapDict[popName] for popName in inPopNames])

#initialize sparse FSs
#first get sets of popus to make FSs for 
if args.FSpops: FSpops = args.FSpops
else:
    FSpops = [[pop] for pop in inPopNames]
    if args.doPairs: FSpops += list(itertools.combinations(inPopNames,2))
    if args.doTrios: FSpops += list(itertools.combinations(inPopNames,3))
    if args.doQuartets: FSpops += list(itertools.combinations(inPopNames,4))

FSs = [SparseFS(len(popGroup), nIntervals) for popGroup in FSpops]

FSrange = list(range(len(FSs)))

###########################################################################################################

linesDone = 0
sitesAnalysed = 0
snpsCounted = 0

### for each line, check if its a scaf we want
for siteData in genoFileReader.siteBySite():
    
    linesDone += 1
    if linesDone % report == 0: sys.stderr.write("\n{} lines done...\n".format(linesDone))
    
    if (include and siteData["scaffold"] not in include) or (exclude and siteData["scaffold"] in exclude): continue
    
    #we will add a count of 1 for this site for each interval it appears in 
    if intervals:
        addValue = intervals.containsPoint(pos=siteData["position"], chrom=siteData["scaffold"])
        #if in no intervals, we move on
        if addValue.sum() == 0: continue
    else: addValue = 1
    
    #check the input data type
    if args.inputType == "genotypes":
        #if it's genotype data, we need to first get the frequencies
        popIndBaseCounts = getPopIndBaseCounts(siteData, args.genoFormat, allSamples, popDict, ploidyDict)
        if popIndBaseCounts is None:
            continue
        
        popBaseCountsArray = getPopBaseCounts(popIndBaseCounts, inPopNames,
                                              subsampleDict=subsampleDict, subsampleIndividuals=args.subsampleIndividuals)
        
        #make sure they all have enough data
        if popBaseCountsArray is None or not np.all(popBaseCountsArray.sum(axis=1) == nHapArray): continue
        
        #if polarizing, get the outgroup base counts
        if outgroup:
            outgroupBaseCounts = popIndBaseCounts[outgroup].sum(axis = 0)
        else: outgroupBaseCounts = None
        
        popTargetCounts = getTargetCounts(popBaseCountsArray, outgroupBaseCounts = outgroupBaseCounts)
        
        if popTargetCounts is None: continue
    
    elif args.inputType == "baseCounts":
        #if it's base counts per population, 
        popBaseCountsArray = np.array([np.array(siteData["GTs"][name].split(","), dtype=float) for name in inPopNames], dtype=int)
        
        #subsample if necessary
        if subsampleDict:
            try: popBaseCountsArray = np.array([downSampleBaseCounts(popBaseCountsArray[i,:], subsampleDict[inPopNames[i]]) for i in range(len(inPopNames))])
            except: continue
        
        if outgroup:
            outgroupBaseCounts = np.array(siteData["GTs"][outgroup].split(","), dtype=float).astype(int)
        else: outgroupBaseCounts = None
        
        popTargetCounts = getTargetCounts(popBaseCountsArray, outgroupBaseCounts = outgroupBaseCounts)
        
        if popTargetCounts is None: continue
        
    else:
        #otherwise assume the input is just counts of the target base
        popTargetCounts = np.array([siteData["GTs"][name] for name in inPopNames])
    
    #if we get here we are going to add this data to the SFS
    popTargetCountsDict = dict(zip(inPopNames, popTargetCounts))

    if args.verbose:
        sys.stderr.write(" ".join(["{}:{}".format(popName,popTargetCountsDict[popName]) for popName in inPopNames]) + "\n")
    
    #for each intervals that this site falls in, make a count
    for i in FSrange:
        #this gets the frequency from each of the relevant pops and uses that list to 'index' the FS
        FSs[i].add([popTargetCountsDict[pop] for pop in FSpops[i]], addValue)
    
    snpsCounted += 1

#write output files. A separate file for each fs
if args.pipe:
    for i in FSrange:
        sys.stdout.write("\n".join(["\t".join([str(x) for x in l]) for l in FSs[i].asChains()]) + "\n")
else:
    for i in FSrange:
        with open(args.pref + "_".join(FSpops[i]) + args.suff, "w") as out:
            out.write("\n".join(["\t".join([str(x) for x in l]) for l in FSs[i].asChains()]) + "\n")
