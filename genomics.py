#functions for manipulating sequences and alignments, working with sliding windows, doing population genetics etx.

import numpy as np
from copy import deepcopy
import sys, string, time, re, math, itertools


##################################################################################################################
#Bits for intyerpreting and manipulating sequence data

DIPLOTYPES = ['A', 'C', 'G', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT', 'AC', 'NN', 'CG', 'AG', 'TT', 'AT', 'CT']
HOMOTYPES = ['A', 'C', 'G', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'N']

diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))
diploHomoDict = dict(zip(DIPLOTYPES,HOMOTYPES))

def haplo(diplo): return diploHaploDict[diplo]

def diplo(pair): return haploDiploDict[pair]

def homo(diplo): return diploHomoDict[diplo]

seqNumDict = {"A":0,"C":1,"G":2,"T":3,"N":-999}

numSeqDict = {0:"A",1:"C",2:"G",3:"T",-999:"N"}


def seqArrayToNumArray(seqArray):
    numArray = np.empty(shape = seqArray.shape, dtype=int)
    for x in ["A","C","G","T","N"]: numArray[seqArray==x] = seqNumDict[x]
    return numArray

def numArrayToSeqArray(numArray):
    seqArray = np.empty(shape = numArray.shape, dtype=str)
    for x in [0,1,2,3,-999]: seqArray[numArray==x] = numSeqDict[x]
    return seqArray


class Genotype:
    def __init__(self, geno, genoFormat, ploidy = None, forcePloidy=False):
        if genoFormat == "phased":
            self.alleles = list(geno)[::2]
            self.phase = geno[1] if len(geno) > 1 and len(geno)%2 == 1 else ""
        elif genoFormat == "alleles" or genoFormat == "pairs" or genoFormat == "haplo":
            self.alleles = list(geno)
            self.phase = ["/"]
        elif genoFormat == "diplo":
            self.alleles = list(haplo(geno))
            self.phase = "/"
        else:
            raise ValueError("Valid genotype formats are 'phased' (eg A/T), 'alleles' (eg AT), 'pairs' (eg AT), 'haplo' (eg A) or 'diplo' (eg W)")
        
        if ploidy is not None:
            ploidyError = ploidy - len(self.alleles)
            if ploidyError != 0:
                if forcePloidy:
                    if ploidyError > 0: self.alleles += ["N"]*ploidyError
                    elif ploidyError < 0:
                        if len(set(self.alleles)) == 1: self.alleles = [self.alleles[0]]*ploidy
                        else: self.alleles = ["N"]*ploidy
                else: raise ValueError("Ploidy doesn't match number of alleles")
        else: ploidy = len(self.alleles)
        
        self.ploidy = ploidy
        
        try: self.numAlleles = np.array([seqNumDict[a] for a in self.alleles])
        except: self.numAlleles = np.array([-999]*self.ploidy)
    
    def isHaploid(self): return self.ploidy == 1
    
    def asPhased(self): return self.phase.join(self.alleles)
    
    def asDiplo(self):
        assert self.ploidy == 2, "Can only convert diploid genotypes to diplotypes."
        return diplo("".join(sorted(self.alleles)))
    
    def asCoded(self, codeDict, missing = None): #code alleles e.g. 0 and 1, with phase (0/1)
        if missing is None: missing = "."
        try: return self.phase.join([codeDict[a] for a in self.alleles])
        except: return self.phase.join([missing]*self.ploidy)
    
    def asCount(self, countAllele, missing = None): # code whole genotype as single value
        if missing is None: missing = -1
        try: return np.bincount(self.numAlleles, minlength = 4)[seqNumDict[countAllele]]
        except: return missing
    
    def isMissing(self): return np.any(self.numAlleles==-999)


#function takes gff file and retrieves coordinates of all CDSs for all mRNAs
def parseGenes(gff):
    #little function to parse the info line
    makeInfoDict = lambda infoString: dict([x.split("=") for x in infoString.strip(";").split(";")])
    output = {}
    for gffLine in gff:
        if len(gffLine) > 1 and gffLine[0] != "#":
            gffObjects = gffLine.split("\t")
            #store all mRNA and CDS data for the particular scaffold
            scaffold = gffObjects[0]
            if scaffold not in output.keys():
                output[scaffold] = {}
            if gffObjects[2] == "mRNA" or gffObjects[2] == "mrna" or gffObjects[2] == "MRNA":
                #we've found a new mRNA
                try: mRNA = makeInfoDict(gffObjects[-1])["ID"]
                except:
                    print gffObjects[-1]
                    raise ValueError("Problem parsing mRNA information.") 
                if mRNA not in output[scaffold].keys():
                    output[scaffold][mRNA] = {'start':int(gffObjects[3]), 'end':int(gffObjects[4]), 'strand':gffObjects[6], 'exons':0, 'cdsStarts':[], 'cdsEnds':[]}
            elif gffObjects[2] == "CDS" or gffObjects[2] == "cds":
                #we're reading CDSs for an existing mRNA
                mRNA = makeInfoDict(gffObjects[-1])["Parent"]
                start = int(gffObjects[3])
                end = int(gffObjects[4])
                output[scaffold][mRNA]['exons'] += 1
                output[scaffold][mRNA]['cdsStarts'].append(start)
                output[scaffold][mRNA]['cdsEnds'].append(end)
    return(output)



#translation table for bases
intab = "ACGTKMRYN"
outtab = "TGCAMKYRN"
trantab = string.maketrans(intab, outtab)


#function to extract a CDS sequence from a genomic sequence given the exon starts, ands and strand
def CDS(seq, seqPos, exonStarts, exonEnds, strand):
  assert len(exonStarts) == len(exonEnds)
  assert len(seq) == len(seqPos)
  seqDict = dict(zip(seqPos,seq))
  cdsSeq = []
  for x in range(len(exonStarts)):
    positions = range(exonStarts[x], exonEnds[x]+1)
    # flip positions if orientation is reverse
    if strand == "-":
      positions.reverse()
      #retrieve bases from seq
    for p in positions:
      try:
        cdsSeq.append(seqDict[p])
      except:
        cdsSeq.append("N")
  cdsSeq = "".join(cdsSeq)  
  #translate if necessary
  if strand == "-":
    cdsSeq = cdsSeq.translate(trantab)
  #if not a multiple of three, remove trailing bases
  overhang = len(cdsSeq) % 3
  if overhang != 0:
    cdsSeq = cdsSeq[:-overhang]
  return cdsSeq


def countStops(cds, includeTerminal=False):
    if includeTerminal:
        triplets = [cds[i:i+3] for i in range(len(cds))[::3]]
    else:
        triplets = [cds[i:i+3] for i in range(len(cds)-3)[::3]]
    stopCount = len([t for t in triplets if t in set(["TAA","TAG","TGA"])])
    return stopCount


#convert one ambiguous sequence into two haploid pseudoPhased sequences
##NOTE this is depricated and should be replaced by splitSeq()
def pseudoPhase(sequence, genoFormat = "diplo"):
    if genoFormat == "pairs": return [[g[0] for g in sequence], [g[1] for g in sequence]]
    elif genoFormat == "phased": return [[g[0] for g in sequence], [g[2] for g in sequence]]
    else:
        pairs = [haplo(g) for g in sequence]
        return [[p[0] for p in pairs], [p[1] for p in pairs]]


def splitSeq(sequence, genoFormat = "phased"):
    assert genoFormat in ["haplo", "diplo", "pairs", "alleles", "phased"]
    if genoFormat == "diplo": sequence = [haplo(d) for d in sequence]
    split = zip(*sequence) 
    #remove phase splitters
    if genoFormat == "phased": split = split[::2]
    return split

#convert a sequence of phased genotypes into two separate sequences
##NOTE this is depricated and should be replaced by splitSeq()
def parsePhase(genotypes):
  first = [geno[0] for geno in genotypes]
  second = [geno[2] for geno in genotypes]
  return [first,second]


#Force diploid sequence to be haploid

def forceHomo(sequence):
    return [homo(s) for s in sequence]



################################################################################################################

#modules for working with individual sites


class GenomeSite:
    
    def __init__(self, genoDict = None, genotypes = None, sampleNames = None, contig = None, position = 0, popDict = {},
                 genoFormat = None, ploidyDict = None, forcePloidy=False):
        #genotypes is a list of genotypes as strings, lists or tuples in any format. e.g. ['AT', 'W', 'T|A', ('A','T)]
        #or use genoDict, which is a dictionary with sample names as the keys. Again, all genotype formats accepted
        if not genoDict:
            assert genotypes is not None, "Either a genotypes dictionary or list must be specified."
            if not sampleNames: sampleNames = map(str, range(len(genotypes)))
            assert len(genotypes) == len(sampleNames), "Genotypes and sample names must be of equal length."
            self.sampleNames = sampleNames
            genoDict = dict(zip(sampleNames, genotypes))
        else:
            self.sampleNames = sorted(genoDict.keys())
        self.contig = contig
        self.position = position
        self.pops = popDict
        self.ploidy = ploidyDict if ploidyDict else dict(zip(sampleNames, [None]*len(sampleNames)))
        
        self.genotypes = {}
        for sample in self.sampleNames:
            self.genotypes[sample] = Genotype(genoDict[sample], genoFormat=genoFormat,
                                              ploidy = self.ploidy[sample],forcePloidy=forcePloidy)
    
    def asList(self, samples = None, pop = None, mode = "phased", alleles = None, codeDict=None, missing=None):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        if mode == "bases":
            return [a for alleles in [self.genotypes[sample].alleles for sample in samples] for a in alleles]
        elif mode == "phased": # like 'A|T' 
            return [self.genotypes[sample].asPhased() for sample in samples]
        elif mode == "diplo": #ACGT and KMRSYW for hets
            return [self.genotypes[sample].asDiplo() for sample in samples]
        elif mode == "alleles": #just bases with no phase
            return [self.genotypes[sample].alleles for sample in samples]
        elif mode == "coded": # vcf format '0/1' - optionally alleles can be provided (REF first)
            if alleles is None: alleles = self.alleles(byFreq = True)
            if codeDict is None: codeDict = dict(zip(alleles, [str(x) for x in range(len(alleles))]))
            return [self.genotypes[sample].asCoded(codeDict, missing) for sample in samples]
        elif mode == "count": # vcf format '0/1' - optionally alleles can be provided (REF first)
            if alleles is None: alleles = self.alleles(byFreq = True)
            countAllele = alleles[-1]
            return [self.genotypes[sample].asCount(countAllele,missing) for sample in samples]
        else:
            raise ValueError("mode must be 'bases', 'phased', 'diplo', 'alleles', 'coded', or 'count'")
    
    def alleles(self, samples = None, pop=None, byFreq = False):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        bases = [a for alleles in [self.genotypes[sample].alleles for sample in samples] for a in alleles]
        alleles, counts = np.unique([b for b in bases if b in "ACGT"], return_counts = True)
        if byFreq: return list(alleles[np.argsort(counts)[::-1]])
        else: return sorted(list(alleles))
    
    def nsamp(self): return len(self.sampleNames)
    
    def changeGeno(self, sample, newGeno, genoFormat="phased"):
        self.genotypes[sample] = Genotype(newGeno, genoFormat=genoFormat,
                                          ploidy = self.ploidy[sample],forcePloidy=forcePloidy)
    
    def hets(self, samples=None):
        if not samples: samples = self.sampleNames
        sampAlleles = self.asList(mode = "alleles")
        sampUniqueAlleles = map(set, sampAlleles)
        nSampAlleles = np.array(map(len, sampUniqueAlleles))
        return 1.*sum(nSampAlleles == 2)/self.nonMissing()
    
    def nonMissing(self, prop=False):
        present = sum([~self.genotypes[sample].isMissing() for sample in self.sampleNames])
        if prop: return 1.*present/(len(sampleNames))
        else: return present
    
    #def plug(self):
    ## plug the major allele in place of missing data
    #popMajor = majorAllele(self.asList(mode="bases"))
    #if len(popMajor) == 1:
        #for sample in popDict[popName]:
            #if site.genotypes[sample].asDiplo() == "N":
                #site.changeGeno(sample, popMajor[0])
    #else:
        #for sample in popDict[popName]:
            #if site.genotypes[sample].asDiplo() == "N":
                #site.changeGeno(sample, random.sample(popMajor,1)[0])


def baseFreqs(bases, asCounts = False, asDict = False):
    counts = np.array([bases.count(i) for i in ["A","C","G","T"]])
    if asCounts: freqs = counts
    else: freqs = counts/sum(counts * 1.)
    if asDict: return dict(zip(["A","C","G","T"], freqs))
    else: return freqs


def majorAllele(bases):
    baseCounts = baseFreqs(bases, asCounts = True, asDict = True)
    m = max(baseCounts.values())
    return [b for b in ["A","C","G","T"] if baseCounts[b] == m]




# method of Wigginton, Cutler and Abecasis, 2005 Am Gen Human Genet. (Adapted from their supplied R code)
def HWEtest(obsHet, obsHom1, obsHom2, side = "both"):
    if obsHom1 < 0 or obsHom2 < 0 or obsHet < 0:
        return -1.0    
    # total genotypes
    N = obsHet + obsHom1 + obsHom2    
    #rare and common number of homozygotes
    obsHomRare,obsHomCom = sorted([obsHom1,obsHom2])    
    #rare allele count
    rare = obsHomRare * 2 + obsHet    
    #initialize probability array
    probs = [0] * (rare + 1)    
    # Find midpoint of the distribution
    mid = math.floor(rare * ( 2 * N - rare) / (2 * N))
    if mid % 2 != rare % 2: mid = mid + 1    
    probs[int(mid)] = 1.0
    mySum = 1.0 
    # Calculate probablities from midpoint down    
    currHet = int(mid)
    currHomRare = int(rare - mid) / 2
    currHomCom = N - currHet - currHomRare    
    while currHet >= 2:
        probs[currHet - 2] = probs[currHet] * currHet * (currHet - 1.0) / (4.0 * (currHomRare + 1.0)    * (currHomCom + 1.0))
        mySum += probs[currHet - 2]        
        # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
        currHet = currHet - 2
        currHomRare = currHomRare + 1
        currHomCom = currHomCom + 1    
    # Calculate probabilities from midpoint up    
    currHet = int(mid)
    currHomRare = int(rare - mid) / 2
    currHomCom = N - currHet - currHomRare    
    while currHet <= rare - 2:
        probs[currHet + 2] = probs[currHet] * 4.0 * currHomRare * currHomCom / ((currHet + 2.0) * (currHet + 1.0))
        mySum += probs[currHet + 2]        
        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        currHet = currHet + 2
        currHomRare = currHomRare - 1
        currHomCom = currHomCom - 1
        
    if side == "top": p = min(1.0, sum(probs[obsHet:(rare+1)]) / mySum)
    elif side == "bottom": p = min(1.0, sum(probs[0:(obsHet+1)]) / mySum)
    else:
        target = probs[obsHet]
        p = min(1.0, sum([prob for prob in probs if prob <= target])/ mySum)
    
    return p


def inHWE(genotypes, P_value, side = "both", verbose = False):
    #genotypes is a list of genotypes as strings, lists or tuples in any format. e.g. ['AT', 'W', 'T|A', ('A','T)]
    site = Site(genotypes = genotypes)
    diplos = site.asList(mode = "diplo")
    diplos = [d for d in diplos if d != "N"]
    if verbose: print diplos
    if len(diplos) == 0: return True
    alleles = site.alleles()
    if len(alleles) == 1: return True
    if len(alleles) > 2: return False
    Hom1Count = int(diplos.count(alleles[0]))
    Hom2Count = int(diplos.count(alleles[1]))
    HetCount = len(diplos) - (Hom1Count + Hom2Count)
    if verbose: print Hom1Count, Hom2Count, HetCount
    p = HWEtest(HetCount,Hom1Count,Hom2Count)
    if verbose: print "P:", p
    if p <= P_value: return False
    else: return True


def siteTest(site,samples=None,minCalls=1,minPopCalls=None,minAlleles=0,maxAlleles=float("inf"),minVarCount=None,maxHet=None,minFreq=None,maxFreq=None,HWE_P=None,HWE_side="both",fixed=False):
    if not samples: samples = site.sampleNames
    #check sufficient number of non-N calls
    if site.nonMissing() < minCalls: return False
    bases = site.asList(mode = "bases", samples=samples)
    baseCounts = baseFreqs(bases, asCounts = True)
    #check min and max alleles 
    nAlleles = len(set(site.alleles(samples)))
    if not minAlleles <= nAlleles <= maxAlleles: return False
    #check variant filters
    if nAlleles > 1:
        # minor allele count
        if minVarCount and sorted(baseCounts)[-2] < minVarCount: return False
        #check maximum heterozygots?
        if maxHet and site.hets(samples) > maxHet: return False
        #if there is a frequency cutoff
        if minFreq and not minFreq <= sorted(baseFreqs(bases))[-2]: return False
        if maxFreq and not sorted(baseFreqs(bases))[-2] <= maxFreq: return False
        #if checking HWE
        if HWE_P:
            #if there are defined pops, check all of them
            if site.pops is not {}:
                for popName in site.pops.keys():
                    if not inHWE(site.asList(pop = popName,mode="diplo"), HWE_P, side = HWE_side): return False
            #otherwise just check all samples
            elif not inHWE(site.asList(mode="diplo"), HWE_P, side = HWE_side): return False
    
    #if there are population-specific filters
    popNames = site.pops.keys()
    if popNames >= 1:        
        for popName in site.pops.keys():
            if minPopCalls:
                popDiplos = [d for d in site.asList(pop=popName, mode = "diplo") if d != "N"]
                if len(popDiplos) < minPopCalls[popName]: return False
    #if we want fixed differences only and there are two or more pops specified
    if fixed:
        #all pops must have only one allele, but taken together must have more than one
        allelesByPop = [site.alleles(pop=popName) for popName in site.pops.keys()]
        if not (set([len(popAlleles) for popAlleles in allelesByPop]) == set([1]) and
                len(set([a for popAlleles in allelesByPop for a in popAlleles])) > 1): return False
    
    #if we get here we've passed all filters
    return True



######################################################################################################################

#modules for working with and analysing alignments

def invertDictOfLists(d):
    new = {}
    for key, lst in d.iteritems():
        for i in lst:
            try: new[i].append(key)
            except: new[i] = [key]
    new
    return new


def makeList(thing):
    if isinstance(thing, basestring): return [thing]
    else:
        try: iter(thing)
        except TypeError: return [thing]
        else: return list(thing)


class Alignment:
    def __init__(self, sequences = None, names=None, groups = None, groupIndDict=None, length = None, numArray = None, positions=None, sampleNames=None):
        assert not sequences is numArray is length is None, "Specify sequences or length of empty sequence object."
        if sequences is not None:
            assert isinstance(sequences, (list,tuple,np.ndarray)), "Sequences must be a list, tuple or numpy array."
            if isinstance(sequences, np.ndarray): seqArray = sequences
            else: seqArray = np.array([list(seq) for seq in sequences])
            
            if numArray is not None: assert numArray.shape == sequences.shape, "Numeric array is different shape from sequence array."
            else: numArray = seqArrayToNumArray(seqArray)
        
        elif numArray is not None:
            assert isinstance(numArray, np.ndarray), "Numeric sequences must be a numpy array."
            seqArray = numArrayToSeqArray(numArray)
        
        else:
            seqArray = np.empty(shape=(0,length), dtype=str)
            numArray = np.empty((0,length), dtype=int)
        
        self.array = seqArray
        self.numArray = numArray
         
        self.nanMask = self.numArray>=0
        
        self.N,self.l = self.array.shape
        
        if positions is not None:
            assert len(positions)==self.l, "Positions must match sequence length."
            self.positions = positions
        else: self.positions = range(1,self.l+1)
        
        if names is None: names = np.arange(self.N)
        else: assert len(names) == self.N, "Incorrect number of names."
        self.names = np.array(names)
        
        if sampleNames is None: sampleNames = self.names
        else: assert len(sampleNames) == self.N, "Incorrect number of sample names."
        self.sampleNames = np.array(sampleNames)
        
        if groups is not None:
            assert len(groups) == self.N, "Incorrect number of groups."
            self.groups = np.array(groups)
            self.indGroupDict = dict(zip(self.names, [makeList(g) for g in self.groups]))
            self.groupIndDict = invertDictOfLists(self.indGroupDict)
        elif groupIndDict is not None:
            self.groupIndDict = groupIndDict
            self.indGroupDict = invertDictOfLists(self.groupIndDict)
            for name in self.names:
                if name not in self.indGroupDict: self.indGroupDict[name] = []
            self.groups = np.array([self.indGroupDict[n] for n in self.names])
        else:
            self.groups = np.array([None]*self.N) #groups is just a list of names, giving the group name for each sample
            self.indGroupDict = dict(zip(self.names, [makeList(g) for g in self.groups]))
            self.groupIndDict = {}
    
    def subset(self, indices = None, names = None, groups = None):
        if indices is None: indices = []
        if names is None: names = []
        if groups is None: groups = []
        names = names + [j for i in [self.groupIndDict[g] for g in groups] for j in i]
        indices += [np.where(self.names == n)[0][0] for n in names]
        indices = np.unique(indices)
        return Alignment(sequences = self.array[indices], numArray=self.numArray[indices],
                        names=self.names[indices], groups=self.groups[indices], sampleNames=self.sampleNames[indices])
    
    def slice(self, indices = None, startPos = None, endPos = None):
        if indices is None:
            if startPos is None: startPos = min(self.positions)
            if endPos is None: endPos = max(self.positions)
            indices = [x for x in range(self.l) if startPos<=self.positions[x]<=endPos]
        
        return Alignment(sequences = self.array[:,indices], names = self.names, groups=self.groups,
                         numArray=self.numArray[:,indices], positions=[self.positions[i] for i in indices])
            
    def column(self,x): return self.array[:,x]
    
    def numColumn(self,x): return self.numArray[:,x]
    
    def distMatrix(self):
        distMat = np.zeros((self.N,self.N))
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                nanMask = self.nanMask[i,:] & self.nanMask[j,:]
                distMat[i,j] = distMat[j,i] = numHamming(self.numArray[i,:][nanMask], self.numArray[j,:][nanMask])
        return distMat
    
    def pairDist(self, i, j):
        nanMask = self.nanMask[i,:] & self.nanMask[j,:]
        return numHamming(self.numArray[i,:][nanMask], self.numArray[j,:][nanMask])
    
    def sampleHet(self, sampleNames=None, asList = False):
        if sampleNames is None: sampleNames,sampleIndices = uniqueIndices(self.sampleNames, preserveOrder=True)
        else: sampleIndices = [np.where(self.sampleNames == sampleName)[0] for sampleName in sampleNames]
        hets = [self.pairDist(x[0],x[1]) if len(x)==2 else np.NaN for x in sampleIndices]
        return dict(zip(sampleNames,hets)) if not asList else hets
    
    def varSites(self, indices=None, names=None):
        if names is not None: indices = np.where(np.in1d(aln.names,names))[0]
        if indices is None: indices = np.arange(self.N)
        return np.where([len(np.unique(self.numArray[indices,x][self.nanMask[indices,x]])) > 1 for x in xrange(self.l)])[0]
    
    def biSites(self): return np.where([len(np.unique(self.numArray[:,x][self.nanMask[:,x]])) == 2 for x in xrange(self.l)])[0]
    
    def siteNonNan(self, sites=None, prop = False):
        if sites is None: sites = range(self.l)
        else: sites = makeList(sites)
        if prop: return np.mean(self.nanMask[:,sites], axis=0)
        return np.sum(self.nanMask[:,sites], axis=0)

    def seqNonNan(self, prop = False):
        if prop: return np.mean(self.nanMask, axis=1)
        return np.sum(self.nanMask, axis=1)
    
    def siteFreqs(self, sites=None, asCounts=False):
        if sites is None: sites = range(self.l)
        else: sites = makeList(sites)
        return np.array([binBaseFreqs(self.numArray[:,x][self.nanMask[:,x]], asCounts=asCounts) for x in sites])
    


def genoToAlignment(seqDict, sampleData=None, genoFormat = "diplo"):
    if sampleData is None: sampleData = SampleData()
    seqNames = []
    sampleNames = []
    groups = []
    haploidSeqs = []
    #first pseudo phase all seqs if necessary
    for indName in seqDict.keys():
        seqList = splitSeq(seqDict[indName], genoFormat)
        ploidy = sampleData.ploidy[indName] if indName in sampleData.ploidy else len(seqList)
        #print seqList
        if ploidy is not None: assert len(seqList) == ploidy, "Sample ploidy (" + str(ploidy) + ") doesn't match number of sequences (" + str(len(seqList)) + ")"
        if ploidy != 1:
            haploidSeqs += seqList
            seqNames += [indName + "_" + letter for letter in string.ascii_uppercase[:ploidy]]
            sampleNames += [indName]*ploidy
            groups += [sampleData.getPop(indName)]*ploidy
        else:
            haploidSeqs.append(forceHomo(seqList[0]))
            seqNames.append(indName)
            sampleNames.append(indName)
            groups.append(sampleData.getPop(indName))
    order = np.argsort(seqNames)
    return Alignment(sequences=[haploidSeqs[i] for i in order],
                     names =   [seqNames[i] for i in order],
                     groups=   [groups[i] for i in order],
                     sampleNames= [sampleNames[i] for i in order])



def binBaseFreqs(numArr, asCounts = False):
    n = len(numArr)
    if n == 0: return np.array([np.NaN]*4)
    else:
        if asCounts: return np.bincount(numArr, minlength=4)
        else: return 1.* np.bincount(numArr, minlength=4) / n


def derivedAllele(inBases, outBases):
    outAlleles = np.unique(outBases)
    inAlleles = np.unique(inBases)
    if len(outAlleles) == 1 and len(inAlleles) == 2 and np.any(outAlleles[0] == inAlleles):
        return inAlleles[inAlleles != outAlleles[0]][0]
    else: return np.nan


def minorAllele(bases):
    alleles = np.unique(bases)
    if len(alleles) == 2:
        alleles, counts = np.unique(bases, return_counts = True)
        return np.random.choice(alleles[counts==min(counts)])
    else: return np.nan


def LD(basesA, basesB, ancA=None, ancB=None):
    arr = np.column_stack([basesA,basesB])
    nanMask = arr >= 0
    goodRows = np.where(np.apply_along_axis(np.all, 1, nanMask))[0]
    arr = arr[goodRows,:]
    N=arr.shape[0]
    allelesA, countsA = np.unique(arr[:,0], return_counts = True)
    allelesB, countsB = np.unique(arr[:,1], return_counts = True)
    if not len(allelesA) == len(allelesB) == 2: return {"D":np.NaN, "Dprime": np.NaN, "r":np.NaN, "r2":np.NaN}
    
    if ancA is None: ancA = allelesA[countsA==max(countsA)][0]
    else: assert ancA in allelesA, "ancestral allele not present"
    
    if ancB is None: ancB = allelesB[countsB==max(countsB)][0]
    else: assert ancB in allelesB, "ancestral allele not present"
    
    boolArr = arr != [ancA,ancB]
    
    pA,pB = np.mean(boolArr, axis = 0)
    pAB = np.mean(np.apply_along_axis(np.all, 1, boolArr))
    
    D = pAB - pA*pB
    Dmin = max([-pA*pB, -(1-pA)*(1-pB)]) if D <0 else min([pA*(1-pB), (1-pA)*pB])
    Dprime = D/Dmin
    r = D / np.sqrt(pA*(1-pA)*pB*(1-pB))
    return {"D":D, "Dprime": Dprime, "r":r, "r2":r**2}


def uniqueIndices(things, preserveOrder = False, asDict=False):
    T,X,I = np.unique(things, return_index=True, return_inverse=True)
    indices = np.array([np.where(I == i)[0] for i in range(len(T))])
    order = np.argsort(X) if preserveOrder else np.arange(len(X))
    return dict(zip(T[order], indices[order])) if asDict else [T[order], indices[order]]

def maxLDphase(aln, sampleIndices=None, stat = "r2"):
    if sampleIndices is None:
        sampleIndices = uniqueIndices(aln.sampleNames, preserveOrder = True)[1]
    
    assert aln.N == sum([len(ind) for ind in sampleIndices]), "Mistmatch between number of indices and sequences"
    assert len(aln.biSites()) == len(aln.varSites()), "Only biallelic or invariant sites are permitted"
    
    nHets = np.array([sum([len(np.unique(aln.numArray[ind,x][aln.nanMask[ind,x]]))>1 for ind in sampleIndices]) for x in range(aln.l)])
    sitesToDo = np.argsort(nHets)[::-1]
    sitesToDo = sitesToDo[nHets[sitesToDo] >= 1]
    
    newNumArray = aln.numArray.copy()
    
    #only do the maximization if there are 2 or more sites with hets
    if len(sitesToDo) >=2:
        #first set the configuration for the first site
        first = sitesToDo[0]
        newNumArray[:,first] = list(itertools.chain(*[sorted(newNumArray[ind, first]) for ind in sampleIndices]))
        newNumArray[:,first] = list(itertools.chain(*[sorted(newNumArray[ind, first]) for ind in sampleIndices]))
        
        for x in range(1,len(sitesToDo)):
            option1 = list(itertools.chain(*[sorted(newNumArray[ind, sitesToDo[x]]) for ind in sampleIndices]))
            option2 = list(itertools.chain(*[sorted(newNumArray[ind, sitesToDo[x]])[::-1] for ind in sampleIndices]))
            LD1 = np.mean([LD(newNumArray[:,sitesToDo[y]],option1)[stat] for y in range(x)])
            LD2 = np.mean([LD(newNumArray[:,sitesToDo[y]],option2)[stat] for y in range(x)])
            newNumArray[:,sitesToDo[x]] = option1 if LD1 >= LD2 else option2
    
    return Alignment(numArray=newNumArray,names=aln.names,groups=aln.groups,
                     positions=aln.positions, sampleNames = aln.sampleNames)



#an older version of sequence distance - using text as opposed to my newer method using numerical arrays
#def seqDistance(seqA, seqB, proportion = True):
    #dist = 0
    #sites = 0
    #for x in xrange(len(seqA)):
        #a,b = seqA[x],seqB[x]
        #if a != "N" and b != "N":
            #sites+=1
            #if a != b:
                #dist += 1
    #if proportion:
        #dist = 1.* dist / sites
    #return dist



## a distance matrix method that uses numerical arrays 
## there is considerable overhead in making the arrays,
## so this isn't fast for pair-wise distance, but is good for sets,
## as you onlyhave to make the array once for each.

def numHamming(numArrayA, numArrayB):
    dif = numArrayA - numArrayB
    return np.mean(dif != 0)


def distMatrix(sequences):
    numSeqs = [[seqNumDict[b] for b in seq] for seq in sequences]
    DNAarray = np.array(numSeqs)
    nanMaskTotal = DNAarray>=0
    N,ln = DNAarray.shape
    distMat = np.zeros((N,N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            nanMask = nanMaskTotal[i,:] & nanMaskTotal[j,:]
            distMat[i,j] = distMat[j,i] = numHamming(DNAarray[i,:][nanMask], DNAarray[j,:][nanMask])
    return distMat


class SampleData:
    def __init__(self, indNames = [], popNames = None, popInds = [], popNumbers = None, ploidyDict = None):
        if popNumbers is None: popNumbers = range(len(popInds))
        if popNames is None: popNames = [str(x) for x in popNumbers]
        assert len(popNames) == len(popInds) == len(popNumbers), "Names, inds and numbers should be same length."
        self.popNames = popNames
        self.popNumbers = popNumbers
        self.popInds = {}
        for x in range(len(popInds)):
            for indName in popInds[x]:
                if indName not in indNames:
                    indNames.append(indName)
            self.popInds[popNames[x]] = popInds[x]
            self.popInds[popNumbers[x]] = popInds[x]
        self.indNames = indNames
        if ploidyDict: self.ploidy = dict(zip(indNames, [ploidyDict[i] for i in indNames]))
        else: self.ploidy = dict(zip(indNames, [2]*len(indNames)))
    
    def getPop(self, indName):
        pop = [p for p in self.popNames if indName in self.popInds[p]]
        if len(pop) == 0: return None
        elif len(pop) == 1: return pop[0]
        else: return tuple(pop)
    
    def getPopNumber(self, popName):
        if popName in self.popNames:
            return self.popNumbers[self.popNames.index(popName)]


def popDiv(Aln, doPairs = True):
    distMat = Aln.distMatrix()
    np.fill_diagonal(distMat, np.NaN) # set all same-with-same to Na
    
    pops,indices = np.unique(Aln.groups, return_inverse = True)
    nPops = len(pops)
    assert nPops > 1, "At least two populations required."
    
    #get population indices - which positions in the alignment correspond to each population
    # this will allow indexing specific pops from the matrix.
    popIndices = [list(np.where(indices==x)[0]) for x in range(nPops)]
    
    output = {}
    
    #pi for each pop
    for x in range(nPops):
        output["pi_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[x])])
    
    if not doPairs: return output

    #pairs
    for x in range(nPops-1):
        for y in range(x+1, nPops):
            #dxy
            output["dxy_" + pops[x] + "_" + pops[y]] = output["dxy_" + pops[y] + "_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[y])])
            
            #fst
            #first get the weightings for each pop
            n_x = len(popIndices[x])
            n_y = len(popIndices[y])
            w = 1.* n_x/(n_x + n_y)
            pi_s = w*(output["pi_" + pops[x]]) + (1-w)*(output["pi_" + pops[y]])
            pi_t = np.nanmean(distMat[np.ix_(popIndices[x]+popIndices[y],popIndices[x]+popIndices[y])])
            output["Fst_" + pops[x] + "_" + pops[y]] = output["Fst_" + pops[y] + "_" + pops[x]] = 1 - pi_s/pi_t
    
    return output

#def ABBABABA(Aln, P1, P2, P3, P4, minData):
    ##subset by population
    #P1Aln = Aln.subset(groups=[P1])
    #P2Aln = Aln.subset(groups=[P2])
    #P3Aln = Aln.subset(groups=[P3])
    #P4Aln = Aln.subset(groups=[P4])
    #P123Aln = Aln.subset(groups=[P1,P2,P3,P4])
    #ABBAsum = BABAsum = maxABBAsum = maxBABAsum = 0.0
    #sitesUsed = 0
    ##get derived frequencies for all biallelic siites
    #for i in P123Aln.biSites():
        #print >> sys.stderr, i
        ##if theres a minimum proportion of sites, check all pops
        #if minData and np.any([A.siteNonNan(i, prop=True) for A in (P1Aln, P2Aln, P3Aln, P4Aln)] < minData): continue
        #allFreqs = Aln.siteFreqs(i)[0] #an array with 4 values, the freq for A,C,G and T
        ## get frequencies for wach pop
        #P1Freqs,P2Freqs,P3Freqs,P4Freqs = [A.siteFreqs(i)[0] for A in (P1Aln, P2Aln, P3Aln, P4Aln)]
        ##check for bad data
        #if np.any(np.isnan(P1Freqs)) or np.any(np.isnan(P2Freqs)) or np.any(np.isnan(P3Freqs)) or np.any(np.isnan(P4Freqs)): continue
        ##if the outgroup is fixed, then that is the ancestral state - otherwise the derived state is the most common allele overall
        #if np.max(P4Freqs) == 1.:
            #anc = np.where(P4Freqs == 1)[0][0] #ancetral allele is which is fixed (get the index)
            #der = [i for i in np.where(allFreqs > 0)[0] if i != anc][0] # derived is the index that is > 0 but not anc
        #else:
            ##der = np.argsort(allFreqs)[-2] # the less common base overall
            #continue
        ##derived allele frequencies
        #P1derFreq = P1Freqs[der]
        #P2derFreq = P2Freqs[der]
        #P3derFreq = P3Freqs[der]
        #P4derFreq = P4Freqs[der]
        #print >> sys.stderr, [P1derFreq, P2derFreq, P3derFreq, P4derFreq]
        #PDderFreq = max(P2derFreq,P3derFreq)
        #print >> sys.stderr, PDderFreq
        ## get weigtings for ABBAs and BABAs
        #ABBA = (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
        #BABA = P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        #maxABBA = (1 - P1derFreq) * PDderFreq * PDderFreq * (1 - P4derFreq)
        #maxBABA = P1derFreq * (1 - PDderFreq) * PDderFreq * (1 - P4derFreq)
        #ABBAsum += (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
        #print >> sys.stderr, "\n"
        #print >> sys.stderr, ABBA
        #print >> sys.stderr, ABBAsum
        #print >> sys.stderr, "\n"
        #BABAsum += P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        #print >> sys.stderr, BABA
        #print >> sys.stderr, BABAsum
        #print >> sys.stderr, "\n"
        #maxABBAsum += (1 - P1derFreq) * PDderFreq * PDderFreq * (1 - P4derFreq)
        #print >> sys.stderr, maxABBA
        #print >> sys.stderr, maxABBAsum
        #print >> sys.stderr, "\n"
        #maxBABAsum += P1derFreq * (1 - PDderFreq) * PDderFreq * (1 - P4derFreq)
        #print >> sys.stderr, maxBABA
        #print >> sys.stderr, maxBABAsum
        #print >> sys.stderr, "\n"
        #sitesUsed += 1
    ##calculate D, fd
    #output = {}
    #try: output["D"] = (ABBAsum - BABAsum) / (ABBAsum + BABAsum)
    #except: output["D"] = np.NaN
    #try:
        #if output["D"] >= 0: output["fd"] = (ABBAsum - BABAsum) / (maxABBAsum - maxBABAsum)
        #else: output["fd"] = np.NaN
    #except: output["fd"] = np.NaN
    #output["ABBA"] = ABBAsum
    #output["BABA"] = BABAsum
    #output["ABBA_P1PDPDO"] = maxABBAsum
    #output["BABA_P1PDPDO"] = maxBABAsum
    #output["sitesUsed"] = sitesUsed
    
    #return output


def ABBABABA(Aln, P1, P2, P3, P4, minData):
    #subset by population
    P1Aln = Aln.subset(groups=[P1])
    P2Aln = Aln.subset(groups=[P2])
    P3Aln = Aln.subset(groups=[P3])
    P4Aln = Aln.subset(groups=[P4])
    P123Aln = Aln.subset(groups=[P1,P2,P3,P4])
    D_numer = D_denom = fd_denom = fdM_denom = 0.0
    sitesUsed = 0
    #get derived frequencies for all biallelic siites
    for i in P123Aln.biSites():
        #if theres a minimum proportion of sites, check all pops
        if minData and np.any([A.siteNonNan(i, prop=True) for A in (P1Aln, P2Aln, P3Aln, P4Aln)] < minData): continue
        allFreqs = Aln.siteFreqs(i)[0] #an array with 4 values, the freq for A,C,G and T
        # get frequencies for wach pop
        P1Freqs,P2Freqs,P3Freqs,P4Freqs = [A.siteFreqs(i)[0] for A in (P1Aln, P2Aln, P3Aln, P4Aln)]
        #check for bad data
        if np.any(np.isnan(P1Freqs)) or np.any(np.isnan(P2Freqs)) or np.any(np.isnan(P3Freqs)) or np.any(np.isnan(P4Freqs)): continue
        #if the outgroup is fixed, then that is the ancestral state - otherwise the derived state is the most common allele overall
        if np.max(P4Freqs) == 1.:
            anc = np.where(P4Freqs == 1)[0][0] #ancetral allele is which is fixed (get the index)
            der = [i for i in np.where(allFreqs > 0)[0] if i != anc][0] # derived is the index that is > 0 but not anc
        else:
            #der = np.argsort(allFreqs)[-2] # the less common base overall
            continue
        #derived allele frequencies
        P1derFreq = P1Freqs[der]
        P2derFreq = P2Freqs[der]
        P3derFreq = P3Freqs[der]
        P4derFreq = P4Freqs[der]
        # get weigtings for ABBAs and BABAs
        D_numer += (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        D_denom += (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq) + P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        
        #fd
        if P2derFreq <= P3derFreq:
            fd_denom += (1 - P1derFreq) * P3derFreq * P3derFreq * (1 - P4derFreq) - P1derFreq * (1 - P3derFreq) * P3derFreq * (1 - P4derFreq)
        else:
            fd_denom += (1 - P1derFreq) * P2derFreq * P2derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P2derFreq * (1 - P4derFreq)
        
        #fdM
        if P2derFreq >= P1derFreq:
            if P2derFreq <= P3derFreq:
                fdM_denom += (1 - P1derFreq) * P3derFreq * P3derFreq * (1 - P4derFreq) - P1derFreq * (1 - P3derFreq) * P3derFreq * (1 - P4derFreq)
            else:
                fdM_denom += (1 - P1derFreq) * P2derFreq * P2derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P2derFreq * (1 - P4derFreq)
        else:
            if P1derFreq <= P3derFreq:
                fdM_denom -= (1 - P3derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq) - P3derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
            else:
                fdM_denom -= (1 - P1derFreq) * P2derFreq * P1derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P1derFreq * (1 - P4derFreq)
        
        sitesUsed += 1
    #calculate D, fd
    output = {}
    try: output["D"] = D_numer / D_denom
    except: output["D"] = np.NaN
    try:
        if output["D"] >= 0: output["fd"] = D_numer / fd_denom
        else: output["fd"] = np.NaN
    except: output["fd"] = np.NaN
    try:
        output["fdM"] = D_numer / fdM_denom
    except: output["fd"] = np.NaN
    output["sitesUsed"] = sitesUsed
    
    return output


##rewriting ABBABABA code to match Hannes Svardal's methods based on numpy arrays

#def ABBABABA(Aln, P1, P2, P3, P4, minData):
    ##subset by population
    #all4Aln = Aln.subset(groups=[P1,P2,P3,P4])
    #P1Aln = all4Aln.subset(groups=[P1])
    #P2Aln = all4Aln.subset(groups=[P2])
    #P3Aln = all4Aln.subset(groups=[P3])
    #P4Aln = all4Aln.subset(groups=[P4])
    ##get derived frequencies for all biallelic siites
    #biSites = all4Aln.biSites()
    #goodSites = (P1Aln.siteNonNan()*1./P1Aln.N >= minData &
                 #P2Aln.siteNonNan()*1./P2Aln.N >= minData &
                 #P3Aln.siteNonNan()*1./P3Aln.N >= minData &
                 #P4Aln.siteNonNan()*1./P4Aln.N >= minData)
    
    #P1freqs = P1Aln.siteFreqs()
    
    


def popSiteFreqs(aln, minData = 0):
    #get population indices
    pops,indices = np.unique(aln.groups, return_inverse = True)
    #subset by population
    popAlns = [aln.subset(groups=pop) for pop in pops]
    #site freqs fro each pop
    _popSiteFreqs = [a.siteFreqs() for a in popAlns]
    #if masking for mising data
    if minData > 0:
        #proportion of inds with non-missing data
        popPropData = [a.siteNonNan for a in popAlns] 
        popDataMask = [propData > minData for propData in popPropData]
        for x in range(len(pops)): _popSiteFreqs[x][~popDataMask[x],:] = np.array([np.nan]*4)
    return _popSiteFreqs


################################################################################################

#modules for working with windows

##Window object class, stores names, sequences and window information


class GenoWindow:
    def __init__(self, scaffold = None, limits=[-np.inf,np.inf],sites = None, names = None, positions = None, ID = None):
        if sites is not None and positions is not None:
            if not len(sites) == len(positions) == 0:
                assert len(set([len(site) for site in sites])) == 1, "Number of genotypes per site must be equal."
                assert len(sites[0]) == len(names), "Number of names must match number of genotypes per site."
                assert len(positions) == len(sites), "Positions must match number of sites"
        else:
            sites = []
            positions = []
        self.names = names if names is not None else []
        self.n = len(self.names)
        self.sites = sites if sites is not None else []
        self.positions = positions if positions is not None else []
        self.scaffold = scaffold
        self.limits = limits
        self.ID = ID
    
    def copy(self): return GenoWindow(scaffold=self.scaffold, limits=self.limits[:],
                                      sites=self.sites[:], names=self.names[:], positions=self.positions[:], ID=self.ID)
    
    ##method for adding
    def addBlock(self, sites, positions):
        assert len(set([len(site) for site in sites])) == 1, "Number of genotypes per site must be equal."
        assert len(sites[0]) == self.n, "Number of genotypes per site must match number of names."
        assert len(positions) == len(sites), "Positions must match number of sites"
        assert np.all(self.limits[0] <= positions <= self.limits[1]), "Position outside of window limit"
        self.sites += sites
        self.positions += positions
    
    def addSite(self, GTs, position=np.NaN, ignorePosition=False):
        assert len(GTs) == self.n, "Number of genotypes per site must match number of names."
        if not ignorePosition:
            assert self.limits[0] <= position <= self.limits[1], "Position: " + str(position) + " outside of window limits: " + "-".join([str(l) for l in self.limits])
        else: position = np.NaN
        self.positions.append(position)
        self.sites.append(GTs)
    
    def seqLen(self): return len(self.positions)
    
    def firstPos(self): return min(self.positions)
    
    def lastPos(self): return max(self.positions)
    
    def slide(self,step=None,newLimits=None):
        #function to slide window along scaffold
        assert step != None or newLimits != None 
        if step: self.limits = [l+step for l in self.limits]
        else: self.limits = newLimits
        i = 0
        while i < len(self.positions) and self.positions[i] < self.limits[0]: i += 1
        #slide positions and sites
        self.positions = self.positions[i:]
        self.sites = self.sites[i:]
    
    def trim(self,right=False,remove=None,leave=None):
        assert remove != None or leave != None
        if not remove: remove=self.seqLen() - leave
        if not right:
            #trim positions and sites
            self.positions = self.positions[remove:]
            self.sites = self.sites[remove:]
        else:
            self.positions = self.positions[:-remove]
            self.sites = self.sites[:-remove]
    
    def seqDict(self, names=None):
        if names is None: names = self.names
        indices = [self.names.index(n) for n in names]
        return dict(zip(names, [[site[i] for site in self.sites] for i in indices]))
    
    def midPos(self):
        try: return int(round(sum(self.positions)/len(self.positions)))
        except: return np.NaN



#site object class for storing the information about a single site
class Site:
    def __init__(self,scaffold=None, position=None, GTs=[]):
        self.scaffold = scaffold
        self.position = position
        self.GTs = GTs

#function to parse a clls line into the Site class
def parseGenoLine(line, splitPhased = False):
    objects = line.split()
    if len(objects) >= 3: site = Site(scaffold = objects[0], position = int(objects[1]), GTs = objects[2:])
    else: site = Site()
    if splitPhased: site.GTs = [a for GT in site.GTs for a in re.split('/|\|', GT)]
    return site


#sliding window generator function
def slidingCoordWindows(genoFile, windSize, stepSize, names = None, splitPhased=False, ploidy = None,
                        include = None, exclude = None, skipDeepcopy = False):
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if names is None: names = allNames
    if splitPhased:
        if ploidy is None: ploidy = [2]*len(allNames)
        ploidyDict = dict(zip(allNames, ploidy))
        #if splitting phased, we need to split names too
        allNames = [n + "_" + letter for n in allNames for letter in string.ascii_uppercase[:ploidyDict[n]]]
        names = [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]
    #indices of samples
    nameIndices = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #window counter
    windowsDone = 0
    #initialise an empty window
    window = GenoWindow()
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    while line:
        #build window
        while site.scaffold == window.scaffold and site.position <= window.limits[1]:
            if site.position >= window.limits[0]:
                #add this site to the window
                window.addSite(GTs=[site.GTs[nameIndices[name]] for name in names], position=site.position)
            #read next line
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        '''if we get here, the line in hand is incompatible with the currrent window
            If the window is not empty, yield it'''
        
        if window.scaffold is not None:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield window.copy()
        
        #now we need to make a new window
        #if on same scaffold, just slide along
        if site.scaffold == window.scaffold:
            window.slide(step = stepSize)
            window.ID = windowsDone + 1
        
        #otherwise we're on a new scaffold (or its the end of the file)
        else:
            #if its one we want to analyse, start new window
            if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
                window = GenoWindow(scaffold = site.scaffold, limits=[1,windSize], names = names, ID = windowsDone + 1)
            
            #if its a scaf we don't want, were going to read lines until we're on one we do want
            else:
                badScaf = site.scaffold
                while site.scaffold == badScaf or (include and site.scaffold not in include and site.scaffold is not None) or (exclude and site.scaffold in exclude and site.scaffold is not None):
                    line = genoFile.readline()
                    site = parseGenoLine(line,splitPhased)
            
        #if we've reached the end of the file, break
        if len(line) <= 1:
            break
    


#sliding window generator function
def slidingSitesWindows(genoFile, windSites, overlap, maxDist = np.inf, minSites = None, names = None,
                        splitPhased=False, ploidy=None, include = None, exclude = None, skipDeepcopy = False):
    if not minSites: minSites = windSites #if minSites < eindSites, windows at ends of scaffolds can still be emmitted
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if names is None: names = allNames
    if splitPhased:
        if ploidy is None: ploidy = [2]*len(allNames)
        ploidyDict = dict(zip(allNames, ploidy))
        #if splitting phased, we need to split names too
        allNames = [n + "_" + letter for n in allNames for letter in string.ascii_uppercase[:ploidyDict[n]]]
        names = [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]
    #indices of samples
    nameIndices = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #window counter
    windowsDone = 0
    #initialise an empty window
    window = GenoWindow()
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    while line:
        #build window
        while site.scaffold == window.scaffold and window.seqLen() < windSites and (window.seqLen() == 0 or site.position - window.firstPos() <= maxDist):
            #add this site to the window
            window.addSite(GTs=[site.GTs[nameIndices[name]] for name in names], position=site.position)
            #read next line
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        '''if we get here, either the window is full, or the line in hand is incompatible with the currrent window
            If the window has more than minSites, yield it'''
        
        
        if window.seqLen() >= minSites:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield deepcopy(window)
            
            #now we need to make a new window
            #if on same scaffold, just trim
            if site.scaffold == window.scaffold:
                window.trim(leave = overlap)
                window.ID = windowsDone + 1
            
            #otherwise we're on a new scaffold (or its the end of the file)
            else:
                #if its one we want to analyse, start new window
                if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
                    window = GenoWindow(scaffold = site.scaffold, names = names, ID = windowsDone + 1)
                
                #if its a scaf we don't want, were going to read lines until we're on one we do want
                else:
                    badScaf = site.scaffold
                    while site.scaffold == badScaf or (include and site.scaffold not in include and site.scaffold is not None) or (exclude and site.scaffold in exclude and site.scaffold is not None):
                    
                        line = genoFile.readline()
                        site = parseGenoLine(line,splitPhased)
        
        #If there are insufficient sites, and we're on the same scaffold, just trim off the furthest left site
        else:
            if site.scaffold == window.scaffold:
                window.trim(remove = 1)
            
            #If we're on a new scaffold, we do as above
            else:
                #if its one we want to analyse, start new window
                if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
                    window = GenoWindow(scaffold = site.scaffold, names = names, ID = windowsDone + 1)
                
                #if its a scaf we don't want, were going to read lines until we're on one we do want
                else:
                    badScaf = site.scaffold
                    while site.scaffold == badScaf or (include and site.scaffold not in include and site.scaffold is not None) or (exclude and site.scaffold in exclude and site.scaffold is not None):
                    
                        line = genoFile.readline()
                        site = parseGenoLine(line,splitPhased)
        
        #if we've reached the end of the file, break
        if len(line) <= 1:
            break


#window generator function using pre-defined coordinates
def predefinedCoordWindows(genoFile, windCoords, names = None, splitPhased=False, ploidy=None, skipDeepcopy = False):
    #get the order of scaffolds
    allScafs = [w[0] for w in windCoords]
    scafs = sorted(set(allScafs), key=lambda x: allScafs.index(x))
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if names is None: names = allNames
    if splitPhased:
        if ploidy is None: ploidy = [2]*len(allNames)
        ploidyDict = dict(zip(allNames, ploidy))
        #if splitting phased, we need to split names too
        allNames = [n + "_" + letter for n in allNames for letter in string.ascii_uppercase[:ploidyDict[n]]]
        names = [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]
    #indices of samples
    nameIndices = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    window = None
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    #we're going to read one line each loop and only releae windows when they're done
    for w in range(len(windCoords)):
        
        #make new window, or if on same scaf just slide it
        ID = windCoords[w][3] if len(windCoords[w]) > 3 else "NA"
        if window and window.scaffold == windCoords[w][0]:
            window.slide(newLimits=[windCoords[w][1],windCoords[w][2]])
            window.ID = ID
        else: window = GenoWindow(scaffold = windCoords[w][0], limits=[windCoords[w][1],windCoords[w][2]], names = names, ID = ID)
        
        #now we need to check that our line in the genome is in a good position
        #if its above the current window - keep reading
        #if it's below the current window - try the next
        
        windScafIdx = scafs.index(window.scaffold)
        
        #if the current scaffold is not in the windows, or is above the window, keep reading
        while site.scaffold and (site.scaffold not in scafs or scafs.index(site.scaffold) < windScafIdx):
            badScaf = site.scaffold
            while site.scaffold == badScaf:
                    line = genoFile.readline()
                    site = parseGenoLine(line,splitPhased)
        
        #if we're on the right scaffold but abve thwe windiow, keep reading
        while site.scaffold == window.scaffold and site.position < window.limits[0]:
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        #if we are in a window - build it
        while site.scaffold == window.scaffold and window.limits[0] <= site.position <= window.limits[1]:
            #add this site to the window
            window.addSite(GTs=[site.GTs[nameIndices[name]] for name in names], position=site.position)
            #read next line
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        '''When we get here, either:
            We're on the right scaffold but below the current window
            We're on a scaffold in the list but below the current window
            We've reached the end of the file
            So we have to yield the current window'''
        
        if skipDeepcopy: yield window
        else: yield deepcopy(window)
        
        if len(line) <= 1: break


#function to read blocks of n lines
#sliding window generator function
def nonOverlappingSitesWindows(genoFile, windSites, names = None, splitPhased=False, ploidy=None, include = None, exclude = None):
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if names is None: names = allNames
    if splitPhased:
        if ploidy is None: ploidy = [2]*len(allNames)
        ploidyDict = dict(zip(allNames, ploidy))
        #if splitting phased, we need to split names too
        allNames = [n + "_" + letter for n in allNames for letter in string.ascii_uppercase[:ploidyDict[n]]]
        names = [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]
    #indices of samples
    nameIndices = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #blocks counter
    windowsDone = 0
    #initialise an empty block
    window = None
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line, splitPhased)
    while True:
        #initialise window
        #if its a scaffold we want to analyse, start new window
        if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
            window = GenoWindow(scaffold = site.scaffold, names = names, ID = windowsDone + 1)
            
        #if its a scaf we don't want, were going to read lines until we're on one we do want
        else:
            window = None
            badScaf = site.scaffold
            while site.scaffold == badScaf or (include and site.scaffold not in include and site.scaffold is not None) or (exclude and site.scaffold in exclude and site.scaffold is not None):
            
                line = genoFile.readline()
                site = parseGenoLine(line,splitPhased)

        #build window
        while window and site.scaffold == window.scaffold and window.seqLen() < windSites:
            #add this site to the window
            window.addSite(GTs=[site.GTs[nameIndices[name]] for name in names], position=site.position)
            #read next line
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        '''if we get here, either the window is full, or the line in hand is incompatible with the currrent window
            If the window has more than minSites, yield it'''
        
        if window:
            windowsDone += 1
            yield window
                
        #if we've reached the end of the file, break
        if len(line) <= 1: break


#function to read entire genoFile into a window-like object
def parseGenoFile(genoFile, names = None, includePositions = False, splitPhased=False, ploidy=None, headerLine = None):
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if names is None: names = allNames
    if splitPhased:
        if ploidy is None: ploidy = [2]*len(allNames)
        ploidyDict = dict(zip(allNames, ploidy))
        #if splitting phased, we need to split names too
        allNames = [n + "_" + letter for n in allNames for letter in string.ascii_uppercase[:ploidyDict[n]]]
        names = [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]
    #indices of samples
    nameIndices = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #initialise an empty window
    window = GenoWindow(names = names)
    for line in iter(genoFile.readline,''):
        site = parseGenoLine(line,splitPhased)
        window.addSite(GTs=[site.GTs[nameIndices[name]] for name in names], position=site.position, ignorePosition= not includePositions)
    
    return window


##########################################################################################################

#functions to make and parse alignment strings in fasta or phylip format

def subset(things,subLen):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]


def makeAlnString(names=None, seqs=None, seqDict=None, outFormat="phylip", lineLen=None):
    assert outFormat=="phylip" or outFormat=="fasta"
    if seqDict: names, seqs = zip(*seqDict.items())
    else: assert len(names) == len(seqs)
    seqs = ["".join(s) for s in seqs]
    output = []
    nSamp = len(names)
    seqLen = max(map(len,seqs))
    if lineLen: seqs = ["\n".join(subset(s,lineLen)) for s in seqs]
    if outFormat == "phylip":
        output.append(" " + str(nSamp) + " " + str(seqLen))
        for x in range(nSamp):
            output.append(names[x] + "   " + seqs[x])
    elif outFormat == "fasta":
        for x in range(nSamp):
            output.append(">" + names[x])
            output.append(seqs[x])
    
    return "\n".join(output) + "\n"


#code to parse alignment strings

def parseFasta(string):
    splitString = string.split(">")[1:]
    names = [s.split()[0] for s in splitString]
    seqs = [s[s.index("\n"):].replace("\n","").replace(" ","") for s in splitString]
    return (names,seqs)


def parsePhylip(string):
    cutString= string[string.index("\n"):]
    names = [l.split()[0] for l in cutString.split("\n") if len(l.split()) == 2]
    for name in names: cutString = cutString.replace("\n" + name + " ", "name")
    splitString = cutString.split("name")[1:]
    seqs = [s.replace("\n","").replace(" ","") for s in splitString]
    return (names,seqs)


############### working with fai

def parseFai(faiFileHandle):
    scafs = []
    lengths = []
    faiLines = fai.readlines()
    for line in faiLines:
        scaf,length,x,y,z = line.split()
        scafs.append(scaf)
        lengths.append(int(length))
    return (tuple(scafs), tuple(lengths),)
