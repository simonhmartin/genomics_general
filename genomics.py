#functions for manipulating sequences and alignments, working with sliding windows, doing population genetics etx.

import numpy as np
import math
from copy import deepcopy
from string import maketrans
import time


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


class Genotype:
    def __init__(self, geno, haploid = False):
        if len(geno) == 1 and geno in DIPLOTYPES:
            self.alleles = list(haplo(geno))
            self.phase = "/"
        elif len(geno) == 2 and "".join(sorted(geno)) in PAIRS:
            self.alleles = list(geno)
            self.phase = "/"
        elif len(geno) == 3 and geno[1] in ["|","/"] and "".join(sorted(geno.split(geno[1]))) in PAIRS:
            self.alleles = geno.split(geno[1])
            self.phase = geno[1]
        else:
            raise ValueError("Unrecognised Genotype")
        
        if haploid:
            assert self.alleles[0] == self.alleles[1], "Biallelic genotype cannot be assigned as haploid"
            self.alleles = self.alleles[:1]
            self.phase = None
            self.ploidy = 1
        else:
            self.ploidy = 2
        
    def isHaploid(self): return self.ploidy == 1
    
    def asPhased(self):
        if self.isHaploid(): return self.alleles
        else: return self.phase.join(self.alleles)
    
    def asDiplo(self):
        if self.isHaploid(): return self.alleles[0]
        else: return diplo("".join(sorted(self.alleles)))




#function to take raw genomic sequence, and the coordinates of exons, and export the concatenated CDS 

#translation table for bases
intab = "ACGTKMRYN"
outtab = "TGCAMKYRN"
trantab = maketrans(intab, outtab)


def CDS(seq, seqPos, exonStarts, exonEnds, strand):
  assert len(exonStarts) == len(exonEnds)
  cdsSeq = []
  for x in range(len(exonStarts)):
    positions = range(exonStarts[x] - 1, exonEnds[x])
    # flip numbers if orientation is reverse
    if strand == "-":
      numbers.reverse()
      #retrieve bases from seq
    for p in positions:
      if p in seqPos:
        cdsSeq.append(seq[seqPos.index[p]])
      else:
        cdsSeq.append("N")
  cdsSeq = "".join(cdsSeq)  
  if strand == "-":
    cdsSeq = cdsSeq.translate(trantab)
  return cdsSeq


def countStops(cds, includeTerminal=False):
    if includeTerminal:
        triplets = [cds[i:i+3] for i in range(len(cds))[::3]]
    else:
        triplets = [cds[i:i+3] for i in range(len(cds)-3)[::3]]
    stopCount = len([t for t in triplets if t in set(["TAA","TAG","TGA"])])
    return stopCount


#convert one ambiguous sequence into two haploid pseudoPhased sequences

def pseudoPhase(sequence, genoFormat = "diplo"):
    if genoFormat == "pairs": pairs = sequence
    else:  pairs = [haplo(s) for s in sequence]
    return [[p[0] for p in pairs], [p[1] for p in pairs]]


#convert a sequence of phased genotypes into two separate sequences

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
    
    def __init__(self, genoDict = None, genotypes = None, sampleNames = None, contig = None, position = 0, popDict = {}, ploidyDict = None):
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
        if ploidyDict:
            self.ploidy = ploidyDict
        else:
            self.ploidy = {}
            for sample in self.sampleNames:
                self.ploidy[sample] = 2
        self.genotypes = {}
        for sample in self.sampleNames:
            self.genotypes[sample] = Genotype(genoDict[sample], haploid = self.ploidy[sample] == 1)
    
    def asList(self, samples = None, pop = None, mode = "phased"):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        if mode == "bases":
            return sum([self.genotypes[sample].alleles for sample in samples],[])
        elif mode == "phased":
            return [self.genotypes[sample].asPhased() for sample in samples]
        elif mode == "diplo":
            return [self.genotypes[sample].asDiplo() for sample in samples]
        elif mode == "alleles":
            return [self.genotypes[sample].alleles for sample in samples]
        else:
            raise ValueError("mode myst be 'bases', 'phased', 'diplo' or 'alleles'")
    
    def alleles(self, samples = None, pop=None):
        if pop: samples = self.pops[pop]
        return sorted(list(set([b for b in self.asList(samples = samples, mode = "bases") if b != "N"])))
    
    def nsamp(self): return len(self.sampleNames)
    
    def changeGeno(self, sample, newGeno):
        self.genotypes[sample] = Genotype(newGeno, haploid = self.ploidy[sample] == 1)
    
    def hets(self, samples=None):
        if not samples: samples = self.sampleNames
        sampAlleles = self.asList(mode = "alleles")
        sampUniqueAlleles = map(set, sampAlleles)
        nSampAlleles = np.array(map(len, sampUniqueAlleles))
        return 1.*sum(nSampAlleles == 2)/self.nNonN()
    
    def nNonN(self):
        return len([d for d in self.asList(mode="diplo") if d != "N"])
    
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


def hets(genotypes):
    #genotypes is a list of genotypes as strings, lists or tuples in any format. e.g. ['AT', 'W', 'T|A', ('A','T)]
    site = Site(genotypes = genotypes)
    sampAlleles = site.asList(mode = "alleles")
    sampUniqueAlleles = map(set, sampAlleles)
    nSampAlleles = np.array(map(len, sampUniqueAlleles))
    return 1.*sum(nSampAlleles == 2)/site.nNonN()


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
    diplos = [d for d in site.asList(mode = "diplo", samples=samples) if d != "N"]
    #check sufficient number of non-N calls
    if len(diplos) < minCalls: return False
    bases = [base for base in site.asList(mode = "bases", samples=samples)]
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
                    if not inHWE(site.asList(pop = popName), HWE_P, side = HWE_side): return False
            #otherwise just check all samples
            elif not inHWE(diplos, HWE_P, side = HWE_side): return False
    
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
        if not set([len(site.alleles(pop=popName)) for popName in site.pops.keys()]) == set([1]) and len(site.alleles(samples = sum([site.pops[popName] for popName in popNames],[]))) > 1 : return False
    
    #if we get here we've passed all filters
    return True



######################################################################################################################

#modules for working with and analysing alignments

numSeqDict = {"A":0,"C":1,"G":2,"T":3,"N":np.NaN}

class Alignment:
    def __init__(self, sequences = None, names=None, groups = None, length = None, numArray = None):
        assert sequences is not None or length is not None, "Specify either sequences or length of empty sequence object."
        if sequences is not None:
            assert isinstance(sequences, (list,tuple,np.ndarray)), "Sequences must be a list, tuple or numpy array."
            if isinstance(sequences, np.ndarray): self.array = sequences
            else: self.array = np.array([list(seq) for seq in sequences])
        else:
            self.array = np.empty((0,length))
            self.numArray = np.empty((0,length))
        
        if numArray is not None:
            assert numArray.shape == sequences.shape, "Numeric array is different shape from sequence array."
            self.numArray = numArray
        else:
            self.numArray = np.array([[numSeqDict[b] for b in seq] for seq in sequences])
         
        self.nanMask = ~np.isnan(self.numArray)
        
        self.N,self.l = self.array.shape
        
        if not names: names = range(self.N)
        else: assert len(names) == self.N, "Incorrect number of names."
        self.names = names
        
        if not groups: groups = [None]*self.N #groups is just a list of names, giving the group name for each sample
        else: assert len(groups) == self.N, "Incorrect number of groups."
        self.groups = groups
                
    def addSeq(self, sequence, name = None, group = None):
        self.array = np.vstack((self.array, np.array(list(sequence))))
        self.numArray = np.vstack((self.numArray, np.array([numSeqDict[b] for b in sequence])))
        self.N += 1
        if not name: name = self.N
        self.names.append(name)
        self.groups.append(group)
        self.nanMask = ~np.isnan(self.numArray)
    
    def subset(self, indices = [], names = [], groups = []):
        indices = [i for i in xrange(self.N) if i in indices or
                   (self.names[i] and self.names[i] in names) or
                   (self.groups[i] and self.groups[i] in groups)]
        
        new = Alignment(sequences = self.array[indices], numArray=self.numArray[indices],
                        names=[self.names[i] for i in indices], groups=[self.groups[i] for i in indices])
        return new
    
    def subsetBySites(self,indices):
        new = Alignment(sequences = self.array[indices], numArray=self.numArray[indices],
                names=[self.names[i] for i in indices], groups=[self.groups[i] for i in indices])
        return new

    
    def column(self,x): return self.array[:,x]
    
    def numColumn(self,x): return self.numArray[:,x]
    
    def distMatrix(self):
        distMat = np.zeros((self.N,self.N))
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                distMat[i,j] = distMat[j,i] = numHamming(self.numArray[i,:], self.numArray[j,:])
        return distMat
    
    def varSites(self): return np.where([np.unique(self.numArray[:,x][self.nanMask[:,x]]) > 1 for x in xrange(self.l)])[0]
    
    def biSites(self): return np.where([len(np.unique(self.numArray[:,x][self.nanMask[:,x]])) == 2 for x in xrange(self.l)])[0]
    
    def siteNonNan(self, sites=None, prop = False):
        if not sites: sites = xrange(self.l)
        if type(sites) is not list: sites = list(sites)
        if prop: return np.array([1.*sum(self.nanMask[:,x])/self.N for x in sites])
        return np.array([sum(self.nanMask[:,x]) for x in sites])
    
    def siteFreqs(self, sites=None):
        if sites is None: sites = range(self.l)
        if type(sites) is not list: sites = list(sites)
        return np.array([binFreqs(self.numArray[:,x][self.nanMask[:,x]].astype(int)) for x in sites])




def genoToAlignment(seqs, sampleData, genoFormat = "diplo"):
    seqNames = []
    groups = []
    pseudoPhasedSeqs = []
    #first pseudo phase all seqs
    for indName in seqs.keys():
        if sampleData.ploidy[indName] == 2:
            pseudoPhasedSeqs += pseudoPhase(seqs[indName], genoFormat)
            seqNames += [indName + "A", indName + "B"]
            groups += [sampleData.getPop(indName)]*2
        else:
            pseudoPhasedSeqs.append(forceHomo(seqs[indName]))
            seqNames.append(indName)
            groups.append(sampleData.getPop(indName))
    order = [seqNames.index(s) for s in sorted(seqNames)]
    return Alignment(sequences=[pseudoPhasedSeqs[i] for i in order],
                     names =   [seqNames[i] for i in order],
                     groups=   [groups[i] for i in order])



def binFreqs(numArr):
    n = len(numArr)
    if n == 0: return np.array([np.NaN]*4)
    else: return 1.* np.bincount(numArr, minlength=4) / n


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
    return np.nanmean(dif[~np.isnan(dif)] != 0)


def distMatrix(sequences):
    numSeqs = [[numSeqDict[b] for b in seq] for seq in seqs]
    DNAarray = np.array(numSeqs)
    N,ln = DNAarray.shape
    distMat = np.zeros((N,N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            distMat[i,j] = distMat[j,i] = numHamming(DNAarray[i,:], DNAarray[j,:])
    return distMat


class SampleData:
    def __init__(self, indNames = [], popNames = None, popInds = [], popNumbers = None, ploidyDict = None):
        if not popNumbers:
            popNumbers = range(len(popInds))
        if not popNames:
            popNames = [str(x) for x in popNumbers]
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


def popDiv(Aln):
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


def ABBABABA(Aln, P1, P2, P3, P4, minData):
    #subset by population
    P1Aln = Aln.subset(groups=[P1])
    P2Aln = Aln.subset(groups=[P2])
    P3Aln = Aln.subset(groups=[P3])
    P4Aln = Aln.subset(groups=[P4])
    P123Aln = Aln.subset(groups=[P1,P2,P3,P4])
    ABBAsum = BABAsum = maxABBAsum = maxBABAsum = 0.0
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
        else:der = np.argsort(allFreqs)[-2] # the less common base overall
        #derived allele frequencies
        P1derFreq = P1Freqs[der]
        P2derFreq = P2Freqs[der]
        P3derFreq = P3Freqs[der]
        P4derFreq = P4Freqs[der]
        PDderFreq = max(P2derFreq,P3derFreq)
        # get weigtings for ABBAs and BABAs
        ABBAsum += (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
        BABAsum += P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        maxABBAsum += (1 - P1derFreq) * PDderFreq * PDderFreq * (1 - P4derFreq)
        maxBABAsum += P1derFreq * (1 - PDderFreq) * PDderFreq * (1 - P4derFreq)
        sitesUsed += 1
    #calculate D, fd
    output = {}
    try: output["D"] = (ABBAsum - BABAsum) / (ABBAsum + BABAsum)
    except: output["D"] = np.NaN
    try:
        if output["D"] >= 0: output["fd"] = (ABBAsum - BABAsum) / (maxABBAsum - maxBABAsum)
        else: output["fd"] = np.NaN
    except: output["fd"] = np.NaN
    output["ABBA"] = ABBAsum
    output["BABA"] = BABAsum
    output["sitesUsed"] = sitesUsed
    
    return output


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

#Window object class, stores names, sequences and window information

class CoordWindow: 
    def __init__(self, scaffold = None, start = None, end = None, seqs = None, names = None, positions = None, ID = None):
        if not names and not seqs:
            names = []
            seqs = []
        elif not names:
            names = [None]*len(seqs)
        elif not seqs:
            seqs = [[] for name in names]
        assert len(names) == len(seqs)
        if not positions:
            positions = []
        if len(seqs) > 0:
            #print len(seqs[0]), len(positions)
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.names = names
        self.positions = positions
        self.seqs = seqs
        self.n = len(self.names)
        self.ID = ID
    
    #method for adding
    def addBlock(self, seqs, positions):
        assert len(seqs) == self.n # ensure correct number seqs is added
        if len(seqs) > 0:
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        for x in range(len(seqs)):
            self.seqs[x] += seqs[x]
        self.positions += positions
    
    def addSite(self, GTs, position):
        assert len(GTs) == self.n # ensure correct number seqs is added
        for x in range(self.n):
            self.seqs[x].append(GTs[x])
        self.positions.append(position)
    
    def seqLen(self):
        return len(self.positions)
    
    def firstPos(self):
        return min(self.positions)
    
    def lastPos(self):
        return max(self.positions)
    
    def slide(self,step=None,newStart=None,newEnd=None):
        #function to slide window along scaffold
        assert step != None or newStart != None 
        if step:
            newStart = self.start + step
            newEnd = self.end + step            
        self.start = newStart
        if newEnd:
            self.end = newEnd
        #find first position beyon newStart
        i = 0
        while i < len(self.positions) and self.positions[i] < newStart:
            i += 1
        #slide positions
        self.positions = self.positions[i:]
        #slide seqs
        self.seqs = [seq[i:] for seq in self.seqs]
    
    def seqDict(self):
        return dict(zip(self.names,self.seqs))
    
    def midPos(self):
        try:
            return int(round(sum(self.positions)/len(self.positions)))
        except:
            return np.NaN


#sites window class - has a fixed number of sites
class SitesWindow: 
    def __init__(self, scaffold = None, seqs = None, names = None, positions = None, ID = None):
        if not names and not seqs:
            names = []
            seqs = []
        elif not names:
            names = [None]*len(seqs)
        elif not seqs:
            seqs = [[] for name in names]
        assert len(names) == len(seqs)
        if not positions:
            positions = []
        if len(seqs) > 0:
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        self.scaffold = scaffold
        self.names = names
        self.positions = positions
        self.seqs = seqs
        self.n = len(self.names)
        self.ID = ID
    
    #method for adding
    def addBlock(self, seqs, positions):
        assert len(seqs) == self.n # ensure correct number seqs is added
        if len(seqs) > 0:
            assert len(seqs[0]) == len(positions) # ensure correct number of positions is given
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs:
                assert type(seq) is list # added sequences must each be a list
        for x in range(len(seqs)):
            self.seqs[x] += seqs[x]
        self.positions += positions
    
    def addSite(self, GTs, position):
        assert len(GTs) == self.n # ensure correct number seqs is added
        for x in range(self.n):
            self.seqs[x].append(GTs[x])
        self.positions.append(position)
    
    def seqLen(self):
        return len(self.positions)
    
    def firstPos(self):
        return min(self.positions)
    
    def lastPos(self):
        return max(self.positions)
    
    def trim(self,right=False,remove=None,leave=None):
        assert remove != None or leave != None
        if not remove: remove=self.seqLen() - leave
        if not right:
            #trim positions
            self.positions = self.positions[remove:]
            #slide seqs
            self.seqs = [seq[remove:] for seq in self.seqs]
        else:
            self.positions = self.positions[:-remove]
            #slide seqs
            self.seqs = [seq[:-remove] for seq in self.seqs]
    
    def seqDict(self):
        return dict(zip(self.names,self.seqs))
    
    def midPos(self):
        try:
            return int(round(sum(self.positions)/len(self.positions)))
        except:
            pass



#sites window class - has a fixed number of sites
class SimpleWindow: 
    def __init__(self, seqs = None, names = None):
        if not names and not seqs:
            names = []
            seqs = []
        elif not names: names = [None]*len(seqs)
        elif not seqs: seqs = [[] for name in names]
        assert len(names) == len(seqs)
        if len(seqs) > 0:
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs: assert type(seq) is list # added sequences must each be a list
        self.names = names
        self.seqs = seqs
        self.n = len(self.names)
    
    #method for adding
    def addBlock(self, seqs):
        assert len(seqs) == self.n # ensure correct number seqs is added
        if len(seqs) > 0:
            assert len(set([len(seq) for seq in seqs])) == 1 #ensure sequences are equal length
            for seq in seqs: assert type(seq) is list # added sequences must each be a list
        for x in range(len(seqs)): self.seqs[x] += seqs[x]
    
    def addSite(self, GTs):
        assert len(GTs) == self.n # ensure correct number seqs is added
        for x in range(self.n): self.seqs[x].append(GTs[x])
    
    def seqLen(self):
        return len(self.seqs[0])
    
    def trim(self,right=False,remove=None,leave=None):
        assert remove != None or leave != None
        if not remove: remove=self.seqLen() - leave
        #trim seqs
        if not right: self.seqs = [seq[remove:] for seq in self.seqs]
        else: self.seqs = [seq[:-remove] for seq in self.seqs]
    
    def seqDict(self): return dict(zip(self.names,self.seqs))
    

#site object class for storing the information about a single site
class Site:
    def __init__(self,scaffold=None, position=None, GTs=[]):
        self.scaffold = scaffold
        self.position = position
        self.GTs = GTs

#function to parse a clls line into the Site class
def parseGenoLine(line, splitPhased = False):
    site = Site()
    objects = line.split()
    if len(objects) >= 3:
        site.scaffold = objects[0]
        site.position = int(objects[1])
        site.GTs = objects[2:]
        if splitPhased: site.GTs = [a for GT in site.GTs for a in [GT[0],GT[-1]]]
    return site

#sliding window generator function
def slidingCoordWindows(genoFile, windSize, stepSize, names = None, splitPhased=False,
                        include = None, exclude = None, skipDeepcopy = False):
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if not names: names = allNames
    if splitPhased:
        #if splitting phased, we need to split names too
        allNames=[name+"_"+x for name in allNames for x in ["A","B"]]
        names=[name+"_"+x for name in names for x in ["A","B"]]
    columns = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #window counter
    windowsDone = 0
    #initialise an empty window
    window = CoordWindow(names = names, ID = 0)
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    while line:
        #build window
        while site.scaffold == window.scaffold and site.position <= window.end:
            #add this site to the window
            window.addSite(GTs=[site.GTs[columns[name]] for name in names], position=site.position)
            #read next line
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        '''if we get here, the line in hand is incompatible with the currrent window
            If the window is not empty, yeild it'''
        
        if window.scaffold is not None:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield deepcopy(window)
        
        #now we need to make a new window
        #if on same scaffold, just slide along
        if site.scaffold == window.scaffold:
            window.slide(step = stepSize)
            window.ID = windowsDone + 1
        
        #otherwise we're on a new scaffold (or its the end of the file)
        else:
            #if its one we want to analyse, start new window
            if (not include and not exclude) or (include and site.scaffold in include) or (exclude and site.scaffold not in exclude):
                window = CoordWindow(scaffold = site.scaffold, start = 1, end = windSize, names = names, ID = windowsDone + 1)
            
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
                        splitPhased=False, include = None, exclude = None, skipDeepcopy = False):
    if not minSites: minSites = windSites #if minSites < eindSites, windows at ends of scaffolds can still be emmitted
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if not names: names = allNames
    if splitPhased:
        #if splitting phased, we need to split names too
        allNames=[name+"_"+x for name in allNames for x in ["A","B"]]
        names=[name+"_"+x for name in names for x in ["A","B"]]
    columns = dict(zip(names, [allNames.index(name) for name in names])) # records site column for each name
    #window counter
    windowsDone = 0
    #initialise an empty window
    window = SitesWindow(names = names, ID = 0)
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    while line:
        #build window
        while site.scaffold == window.scaffold and window.seqLen() < windSites and (window.seqLen() == 0 or site.position - window.firstPos() <= maxDist):
            #add this site to the window
            window.addSite(GTs=[site.GTs[columns[name]] for name in names], position=site.position)
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
                    window = SitesWindow(scaffold = site.scaffold, names = names, ID = windowsDone + 1)
                
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
                    window = SitesWindow(scaffold = site.scaffold, names = names, ID = windowsDone + 1)
                
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
def predefinedCoordWindows(genoFile, windCoords, names = None, splitPhased=False, skipDeepcopy = False):
    #get the order of scaffolds
    allScafs = [w[0] for w in windCoords]
    scafs = sorted(set(allScafs), key=lambda x: allScafs.index(x))
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if not names: names = allNames
    if splitPhased:
        #if splitting phased, we need to split names too
        allNames=[name+"_"+x for name in allNames for x in ["A","B"]]
        names=[name+"_"+x for name in names for x in ["A","B"]]
    columns = dict(zip(names, [allNames.index(name) for name in names])) # records file column for each name
    #window counter
    w = 0
    window = None
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    #we're going to read one line each loop and only releae windows when they're done
    for w in range(len(windCoords)):
        
        #make new window, or if on same scaf just slide it
        if window and window.scaffold == windCoords[w][0]:
            window.slide(newStart = windCoords[w][1], newEnd = windCoords[w][2])
            window.ID = w
        else: window = CoordWindow(scaffold = windCoords[w][0], start=windCoords[w][1], end=windCoords[w][2], names = names, ID = w)
        
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
        while site.scaffold == window.scaffold and site.position < window.start:
            line = genoFile.readline()
            site = parseGenoLine(line,splitPhased)
        
        #if we are in a window - build it
        while site.scaffold == window.scaffold and window.start <= site.position <= window.end:
            #add this site to the window
            window.addSite(GTs=[site.GTs[columns[name]] for name in names], position=site.position)
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


#function to read entire genoFile into a window-like object
def parseGenoFile(genoFile, names = None,splitPhased=False):
    #get file headers
    headers = genoFile.readline().split()
    allNames = headers[2:]
    if not names: names = allNames
    if splitPhased:
        #if splitting phased, we need to split names too
        allNames=[name+"_"+x for name in allNames for x in ["A","B"]]
        names=[name+"_"+x for name in names for x in ["A","B"]]
    columns = dict(zip(names, [allNames.index(name) for name in names])) # records site column for each name
    #initialise an empty window
    window = SimpleWindow(names = names)
    #read first line
    line = genoFile.readline()
    site = parseGenoLine(line,splitPhased)
    while line:
        window.addSite(GTs=[site.GTs[columns[name]] for name in names])
        #read next line
        line = genoFile.readline()
        site = parseGenoLine(line,splitPhased)
        #if we've reached the end of the file, break
        if len(line) <= 1:
            break
    
    return window


##########################################################################################################

#functions to make and parse alignment strings in fasta or phylip format

def subset(things,subLen):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]


def makeAlnString(names, seqs, format="phylip", lineLen=None):
    assert format=="phylip" or format=="fasta"
    assert len(names) == len(seqs)
    seqs = ["".join(s) for s in seqs]
    output = []
    nSamp = len(names)
    seqLen = max(map(len,seqs))
    if lineLen: seqs = ["\n".join(subset(s,lineLen)) for s in seqs]
    if format == "phylip":
        output.append(" " + str(nSamp) + " " + str(seqLen))
        for x in range(nSamp):
            output.append(names[x] + "   " + seqs[x])
    elif format == "fasta":
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
