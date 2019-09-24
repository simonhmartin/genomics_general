#functions for manipulating sequences and alignments, working with sliding windows, doing population genetics etx.

import numpy as np
from copy import copy, deepcopy
from collections import defaultdict
import sys, string, time, re, math, itertools, random

np.seterr(divide='ignore')

##################################################################################################################

#Bits for intyerpreting and manipulating sequence data

DIPLOTYPES = ('A',  'C',  'G',  'K',  'M',  'N',  'S',  'R',  'T',  'W',  'Y')
PAIRS =      ('AA', 'CC', 'GG', 'GT', 'AC', 'NN', 'CG', 'AG', 'TT', 'AT', 'CT')
HOMOTYPES =  ('A',  'C',  'G',  'N',  'N',  'N',  'N',  'N',  'T',  'N',  'N')

IUPAC =      ('A',  'C',  'G',  'T',  'M',  'R',  'W',  'S',  'Y',  'K',  'V',   'H',   'D',   'B',   'N')
ALLTYPES =   ('A',  'C',  'G',  'T',  'AC', 'AG', 'AT', 'CG', 'CT', 'GT', 'ACG', 'ACT', 'AGT', 'CGT', 'ACGT')

diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))
diploHomoDict = dict(zip(DIPLOTYPES,HOMOTYPES))
basesIupacDict = dict(zip(ALLTYPES,IUPAC))
iupacBasesDict = dict(zip(IUPAC,ALLTYPES))

def haplo(diplo): return diploHaploDict[diplo]

def diplo(pair): return haploDiploDict[pair]

def homo(diplo): return diploHomoDict[diplo]

seqNumDict = {"A":0,"C":1,"G":2,"T":3,"N":-999}

numSeqDict = {0:"A",1:"C",2:"G",3:"T",-999:"N"}


#translation tables - method epends on version

if sys.version_info>=(3,0):
    #translation for conversion of missing bases to gaps
    missingtrans = str.maketrans("Nn", "--")
    #translation table for bases
    complementTrans = str.maketrans("ACGTKMRYVHBDN", "TGCAMKYRBDVHN")
else:
    #translation for conversion of missing bases to gaps
    missingtrans = string.maketrans("Nn", "--")
    #translation table for bases
    complementTrans = string.maketrans("ACGTKMRYVHBDN", "TGCAMKYRBDVHN")

complementDict = dict(zip(list("ACGTKMRYVHBDN"), list("TGCAMKYRBDVHN")))

def complement(seq):
    if type(seq) == str: return seq.translate(complementTrans)
    else: return [complementDict[a] for a in seq]

def revComplement(seq):
    if type(seq) == str: return seq.translate(complementTrans)[::-1]
    else: return [complementDict[a] for a in seq[::-1]]

def allPossibleSeqs(seq, ignoreNs = True):
    if ignoreNs: basesList = [iupacBasesDict[s] if s != "N" else "N" for s in seq]
    else: basesList = [iupacBasesDict[s] for s in seq]
    seqs = [[]]
    for bases in basesList:
        for x in range(len(seqs)):
            seqs[x].append(bases[0])
            for b in bases[1:]:
                newSeq = seqs[x][:]
                newSeq[-1] = b
                seqs.append(newSeq)
    return ["".join(s) for s in seqs]

def seqArrayToNumArray(seqArray):
    numArray = np.empty(shape = seqArray.shape, dtype=int)
    for x in ["A","C","G","T","N"]: numArray[seqArray==x] = seqNumDict[x]
    return numArray

def numArrayToSeqArray(numArray):
    seqArray = np.empty(shape = numArray.shape, dtype=str)
    for x in [0,1,2,3,-999]: seqArray[numArray==x] = numSeqDict[x]
    return seqArray

def alleles(bases):
    s = set(bases)
    return [i for i in "ACGT" if i in s]

################################################################################

#working with coding sequences

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def translate(sequence):
    """Return the translated protein from 'sequence' assuming +1 reading frame"""
    return ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])

def possibleCodons(pos1Alleles,pos2Alleles,pos3Alleles):
    return ["".join(x) for x in itertools.product(pos1Alleles,pos2Alleles,pos3Alleles)]

def possibleAAs(pos1Alleles,pos2Alleles,pos3Alleles):
    return sorted(set([translate(codon) for codon in possibleCodons(pos1Alleles,pos2Alleles,pos3Alleles)]))

#function that takes three sets of allels, corresponding to the alleles at the three codon positions
#and outputs whether each one is synonymous or nonsynonymous
#rejects cases where two sites are variable
def synNon(pos1Alleles,pos2Alleles,pos3Alleles):
    output = ["NA","NA","NA"]
    #first check that each position has at least one allele and at most one position is variable
    nAlleles = [len(alleles) for alleles in (pos1Alleles,pos2Alleles,pos3Alleles,)]
    if not sorted(nAlleles) == [1,1,2]:
        return output
    else:
        focal = nAlleles.index(2)
        output[focal] = "syn" if len(possibleAAs(pos1Alleles,pos2Alleles,pos3Alleles)) == 1 else "non"
    return output


#dictionary to tell you how degenerate a site is based on how many unique amino acids are formed when that site is mutated
#eg if four distinct amino acids can be formed, the site is 0-fold degenerate
degenDict = {4:0, 3:2, 2:2, 1:4}

def degeneracy(pos1Alleles,pos2Alleles,pos3Alleles):
    #if its invariant, then they all get a degeneracy
    if len(pos1Alleles) == len(pos2Alleles) == len(pos3Alleles) == 1:
        output = [degenDict[len(possibleAAs("ACGT",pos2Alleles,pos3Alleles))],
                  degenDict[len(possibleAAs(pos1Alleles,"ACGT",pos3Alleles))],
                  degenDict[len(possibleAAs(pos1Alleles,pos2Alleles,"ACGT"))]]
    
    elif len(pos1Alleles) == 2 and len(pos2Alleles) == len(pos3Alleles) == 1:
        output = [degenDict[len(possibleAAs("ACGT",pos2Alleles,pos3Alleles))], "NA", "NA"]
    
    elif len(pos2Alleles) == 2 and len(pos1Alleles) == len(pos3Alleles) == 1:
        output = ["NA", degenDict[len(possibleAAs(pos1Alleles,"ACGT",pos3Alleles))], "NA"]
    
    elif len(pos3Alleles) == 2 and len(pos1Alleles) == len(pos2Alleles) == 1:
        output = ["NA","NA", degenDict[len(possibleAAs(pos1Alleles,pos2Alleles,"ACGT"))]]
    
    else:
        output = ["NA","NA","NA"]
    
    return output


#function takes gff file and retrieves coordinates of all CDSs for all mRNAs
def parseGenes(gff):
    #little function to parse the info line
    makeInfoDict = lambda infoString: dict([x.split("=") for x in infoString.strip(";").split(";")])
    output = {}
    for gffLine in gff:
        if len(gffLine) > 1 and gffLine[0] != "#":
            gffObjects = gffLine.strip().split("\t")
            #store all mRNA and CDS data for the particular scaffold
            scaffold = gffObjects[0]
            if scaffold not in output.keys():
                output[scaffold] = {}
            if gffObjects[2] == "mRNA" or gffObjects[2] == "mrna" or gffObjects[2] == "MRNA":
                #we've found a new mRNA
                try: mRNA = makeInfoDict(gffObjects[-1])["ID"]
                except:
                    raise ValueError("Problem parsing mRNA information: " + gffObjects[-1]) 
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


#function to extract a CDS sequence from a genomic sequence given the exon starts, ands and strand
def CDSpositions(exonStarts, exonEnds, strand, trim=False):
    nExons = len(exonStarts)
    assert nExons == len(exonEnds)
    
    #for each exon, extract a list of scaffold positions
    codingPositions = [list(range(exonStarts[i], exonEnds[i] + 1)) for i in range(nExons)]
    
    #reverse the positions if necessary
    if strand == "-":
        for _codingPositions_ in codingPositions: _codingPositions_.reverse()
    
    #unlist them to make one long set of positions
    codingPositions = [i for _codingPositions_ in codingPositions for i in _codingPositions_]
    
    if trim:
        overhang = len(codingPositions) % 3
        if overhang != 0: codingPositions = codingPositions[:-overhang]
    
    return codingPositions

#function to extract a CDS sequence from a genomic sequence given the exon starts, ands and strand
def CDSsequence(exonStarts, exonEnds, strand, seqDict=None, seq=None, seqPos=None, trim=True):
    
    if seqDict is None:
        #if dictionary of bases for each position is not provided, make it
        assert len(seq) == len(seqPos)
        seqDict = defaultdict(lambda: "N", zip(seqPos, seq))
    
    #for each exon, extract a list of scaffold positions
    codingPositions = CDSpositions(exonStarts, exonEnds, strand, trim=trim)
    
    cdsSeq = "".join([seqDict[p] for p in codingPositions])
    
    if strand=="-": cdsSeq = cdsSeq.translate(complementTrans)
    
    return cdsSeq


def countStops(cds, includeTerminal=False):
    if includeTerminal:
        triplets = [cds[i:i+3] for i in range(len(cds))[::3]]
    else:
        triplets = [cds[i:i+3] for i in range(len(cds)-3)[::3]]
    stopCount = len([t for t in triplets if t in set(["TAA","TAG","TGA"])])
    return stopCount



################################################################################
########### some general list manipulation code

#subset list into smaller lists
def subset(things,subLen):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]

#similar to above, but can have variable sizes of smaller lists, or can specify the number of chunks
def chunkList(l, nChunks = None, chunkSize = None, return_indices=False):
    N = len(l)
    assert not nChunks is chunkSize is None
    if nChunks is not None:
        assert N % nChunks == 0, "list must be divizable by number of chunks"
        chunkSize = [N/nChunks]*nChunks
    elif isinstance(chunkSize, int):
        assert N % chunkSize == 0, "list must be divizable by chunk size"
        chunkSize = [chunkSize]*(N/chunkSize)
    elif len(chunkSize) == 1:
        assert N % chunkSize[0] == 0, "list must be divizable by chunk size"
        chunkSize*=(N/chunkSize[0])
    else: assert N == sum(chunkSize), "Chunk sizes must sum to list length"
        
    indices = []
    r = range(N)
    i = 0
    for c in chunkSize:
        indices.append(range(i,i+c))
        i = i+c
    if return_indices: return ([[l[x] for x in ind] for ind in indices], indices,)
    else: return [[l[x] for x in ind] for ind in indices]


def invertDictOfLists(d):
    new = {}
    for key, lst in d.items():
        for i in lst:
            try: new[i].append(key)
            except: new[i] = [key]
    new
    return new


def makeList(thing):
    if isinstance(thing, str): return [thing]
    else:
        try: iter(thing)
        except TypeError: return [thing]
        else: return list(thing)

def uniqueIndices(things, preserveOrder = False, asDict=False):
    T,X,I = np.unique(things, return_index=True, return_inverse=True)
    indices = np.array([np.where(I == i)[0] for i in range(len(T))])
    order = np.argsort(X) if preserveOrder else np.arange(len(X))
    return dict(zip(T[order], indices[order])) if asDict else [T[order], indices[order]]

#################################################################################################


class Genotype:
    __slots__ = ['geno', 'genoFormat', 'ploidy', 'forcePloidy', 'alleles', 'phase', 'numAlleles']
    
    def __init__(self, geno, genoFormat, ploidy = None, forcePloidy=False):
        if genoFormat == "phased":
            self.alleles = list(geno)[::2]
            self.phase = geno[1] if len(geno) > 1 and len(geno)%2 == 1 else "/"
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
        
        #now make the alleles immutable
        self.alleles = tuple(self.alleles)
        
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
    
    def asBaseCounts(self):
        return np.bincount(self.numAlleles[self.numAlleles >= 0], minlength = 4)
    
    def isMissing(self): return np.any(self.numAlleles==-999)


#convert one ambiguous sequence into two haploid pseudoPhased sequences
##NOTE this is depricated and should be replaced by splitSeq()
def pseudoPhase(sequence, genoFormat = "diplo"):
    if genoFormat == "pairs": return [[g[0] for g in sequence], [g[1] for g in sequence]]
    elif genoFormat == "phased": return [[g[0] for g in sequence], [g[2] for g in sequence]]
    else:
        pairs = [haplo(g) for g in sequence]
        return [[p[0] for p in pairs], [p[1] for p in pairs]]

def splitSeq(sequence, genoFormat = "phased"):
    assert genoFormat in ("haplo", "diplo", "pairs", "alleles", "phased",)
    if genoFormat == "diplo": sequence = [haplo(d) for d in sequence]
    split = list(zip(*sequence)) 
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


#make haploid sequences N-ploid
def haploToPhased(seqs, seqNames=None, ploidy=2, randomPhase=False):
    #phase if ploidy is not 1
    if ploidy is not 1:
        _ploidy_ = makeList(ploidy)
        Nseqs = len(seqs)
        if len(_ploidy_) == 1:
            assert Nseqs % _ploidy_[0] == 0, "Sequence number must be divizable by ploidy"
            _ploidy_ = _ploidy_*(Nseqs/_ploidy_[0])
        else:
            assert Nseqs == sum(_ploidy_), "Ploidys must sum to number of sequences"
        
        indices = chunkList(range(Nseqs), chunkSize=_ploidy_, return_indices=True)[1]
        
        zipSeqs = [zip(*[seqs[x] for x in ind]) for ind in indices]
        #randomize phase if necessary
        if randomPhase:
            for i in range(len(indices)):
                if _ploidy_[i] > 1:
                    for j in range(len(zipSeqs[i])):
                        zipSeqs[i][j] = random.sample(zipSeqs[i][j], _ploidy_[i])
        seqs = [["|".join(x) for x in zipSeq] for zipSeq in zipSeqs]
        #if seqNames provided, return a tuple
        if seqNames != None:
            assert len(seqNames) == Nseqs, "incorrect number of sequence names"
            seqNames = ["_".join([seqNames[x] for x in ind]) for ind in indices]
            return (seqs, seqNames,)
        #otherwise just return the seqs
        return seqs

def makeHaploidNames(names,ploidy=2):
    ploidy=makeList(ploidy)
    if len(ploidy) == 1: ploidy = ploidy*len(names)
    ploidyDict = dict(zip(names, ploidy))
    return [n + "_" + letter for n in names for letter in string.ascii_uppercase[:ploidyDict[n]]]

def makePhasedNames(names,ploidy=2):
    nameGroups = chunkList(names, chunkSize=ploidy)
    return ["_".join(group) for group in nameGroups]


################################################################################################################

#modules for working with individual sites


class GenomeSite:
    
    def __init__(self, genoDict = None, genotypes = None, sampleNames = None, contig = None, position = 0, popDict = {},
                 genoFormat = None, ploidyDict = None, forcePloidy=False, precompGTs=None, addToPrecomp=True):
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
        self.ploidy = ploidyDict if ploidyDict else dict(zip(self.sampleNames, [None]*len(self.sampleNames)))
        
        self.genotypes = {}
        for sample in self.sampleNames:
            #now, for each sample, if using precomputed genotypes, we check if the geno has already been computed
            if precompGTs and genoDict[sample] in precompGTs[sample]:
                self.genotypes[sample] = precompGTs[sample][genoDict[sample]]
            
            #if not, compute gentype normally
            else:
                self.genotypes[sample] = Genotype(genoDict[sample], genoFormat=genoFormat,
                                                  ploidy = self.ploidy[sample],forcePloidy=forcePloidy)
                
                #add to precomputed
                if precompGTs and addToPrecomp:
                    precompGTs[sample][genoDict[sample]] = self.genotypes[sample]
    
    
    def asList(self, samples = None, pop = None, mode = "phased", alleles = None,
               codeDict=None, missing=None, alleleOrder=None, countAllele=None):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        if mode == "bases":
            #if we want the bases returned in order of their overall frequency
            if alleleOrder == "freq":
                siteAlleles = self.alleles(samples=samples,byFreq = True) + ["N"]
                return [a for sample in samples for a in sorted(self.genotypes[sample].alleles, key=lambda x: siteAlleles.index(x))]
            #otherwise just return the bases as they appear
            else:
                return [a for alleles in [self.genotypes[sample].alleles for sample in samples] for a in alleles]
        elif mode == "alleles": #just bases with no phase
            #if we want the bases returned in order of their overall frequency
            if alleleOrder == "freq":
                siteAlleles = self.alleles(samples=samples,byFreq = True) + ["N"]
                return ["".join(sorted(self.genotypes[sample].alleles,
                                                                 key=lambda x: siteAlleles.index(x))) for sample in samples]
            #otherwise just return the bases as they appear
            else:
                return [self.genotypes[sample].alleles for sample in samples]
        if mode == "numeric":
            return np.concatenate([self.genotypes[sample].numAlleles for sample in samples])
        elif mode == "numAlleles": #numpy array of numeric alleles
            return [self.genotypes[sample].numAlleles for sample in samples]
        elif mode == "phased": # like 'A|T' 
            return [self.genotypes[sample].asPhased() for sample in samples]
        elif mode == "diplo": #ACGT and KMRSYW for hets
            return [self.genotypes[sample].asDiplo() for sample in samples]
        elif mode == "coded": # vcf format '0/1' - optionally alleles can be provided (REF first)
            if alleles is None: alleles = self.alleles(byFreq = True)
            if codeDict is None: codeDict = dict(zip(alleles, [str(x) for x in range(len(alleles))]))
            return [self.genotypes[sample].asCoded(codeDict, missing) for sample in samples]
        elif mode == "count":
            if countAllele is None:
                if alleles is None: alleles = self.alleles(byFreq = True)
                countAllele = alleles[-1]
            return [self.genotypes[sample].asCount(countAllele,missing) for sample in samples]
        else:
            raise ValueError("mode must be 'bases', 'alleles', 'numeric', 'numAlleles', 'phased', 'diplo', 'coded', or 'count'")
    
    def baseFreqs(self, samples = None, pop=None, asCounts=False):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        numBases = np.concatenate([self.genotypes[sample].numAlleles for sample in samples])
        return binBaseFreqs(numBases[numBases >= 0], asCounts = asCounts)
    
    def alleles(self, samples = None, pop=None, byFreq = False, numeric=False):
        if pop: samples = self.pops[pop]
        if not samples: samples = self.sampleNames
        counts = self.baseFreqs(samples=samples,asCounts = True)
        idx = counts>0
        alleles = np.array(["A","C","G","T"])[idx] if not numeric else idx
        counts = counts[idx]
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


def binBaseFreqs(numArr, asCounts = False):
    n = len(numArr)
    if n == 0:
        if asCounts: return np.zeros(4, dtype=int)
        else: return np.array([np.NaN]*4)
    else:
        if asCounts: return np.bincount(numArr, minlength=4)
        else: return 1.* np.bincount(numArr, minlength=4) / n


#quickest method I could find to determine if a numeric array is variable or not
def numVar(numArray):
    return max(np.bincount(numArray)) != len(numArray)

#timeit.timeit("numVar(numArray)", setup="from __main__ import numVar, numArray", number = 1000000)

#site-wise pi for multi-allelic sites (input is counts of four bases)
def baseCountPi(baseCounts):
    N = sum(baseCounts)
    return (baseCounts[0]*baseCounts[1] +
            baseCounts[0]*baseCounts[2] + 
            baseCounts[0]*baseCounts[3] + 
            baseCounts[1]*baseCounts[2] +
            baseCounts[1]*baseCounts[3] + 
            baseCounts[2]*baseCounts[3]) / (.5*N*(N-1))


def TajimaD(n,S,theta_pi):
    # https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
    a= sum( 1./i for i in range(1, n))
    theta_w = 1.*S/a # M 
    a2 = sum( 1./(i**2) for i in range(1, n))
    b1 = (n + 1.) / (3*(n-1))
    b2 = (2. * (n**2 + n + 3)) / (9*n*(n-1))
    c1 = b1 - (1./a)
    c2 = b2 - ((n+2)/(a*n)) + a2/(a**2)
    e1 = c1/a
    e2 = c2/(a**2 + a2)        
    d = theta_pi - theta_w
    D = d / np.sqrt(e1*S + e2*S*(S-1))
    return D



def derivedAllele(inBases=None, outBases=None,
                  inBaseCounts=None, outBaseCounts=None,
                  inAlleles=[], outAlleles=[],
                  maxOneDerivedAllele=True, numeric=False):
    
    if inBases is not None: inAlleles = np.unique(inBases)
    elif inBaseCounts is not None:
        if not numeric: inAlleles=["ACGT"[i] for i in range(4) if inBaseCounts[i] > 0]
        else: inAlleles= np.where(np.array(inBaseCounts) > 0)[0]
    
    if outBases is not None: outAlleles = np.unique(outBases)
    elif outBaseCounts is not None:
        if not numeric: outAlleles=["ACGT"[i] for i in range(4) if outBaseCounts[i] > 0]
        outAlleles= np.where(np.array(outBaseCounts) > 0)[0]
    
    if not isinstance(inAlleles, np.ndarray): inAlleles = np.array(inAlleles)
    if not isinstance(outAlleles, np.ndarray): outAlleles = np.array(outAlleles)
    
    if maxOneDerivedAllele and len(outAlleles) == 1 and len(inAlleles) == 2 and np.any(outAlleles[0] == inAlleles):
        return inAlleles[inAlleles != outAlleles[0]][0]
    
    elif not maxOneDerivedAllele and len(outAlleles) == 1 and len(inAlleles) >= 2 and np.any(outAlleles[0] == inAlleles):
        return inAlleles[inAlleles != outAlleles[0]]
    
    elif numeric or isinstance(inAlleles[0], np.int) or isinstance(inAlleles[0], np.float): return np.nan
    else: return "N"

def minorAllele(bases):
    alleles = np.unique(bases)
    if len(alleles) == 2:
        alleles, counts = np.unique(bases, return_counts = True)
        return np.random.choice(alleles[counts==min(counts)])
    else: return np.nan


def consensus(bases):
    x = "".join(np.unique([b for b in bases if b in "ACGT"]))
    if x == "": x = "ACGT"
    return(basesIupacDict[x])


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

def inHWE(diplos, P_value, side = "both", verbose = False):
    diplos = [d for d in diplos if d != "N"]
    if verbose: sys.stderr.write(diplos)
    if len(diplos) == 0: return True
    alleles = unique([haplo(d) for d in diplos])
    if len(alleles) == 1: return True
    if len(alleles) > 2: return False
    Hom1Count = int(diplos.count(alleles[0]))
    Hom2Count = int(diplos.count(alleles[1]))
    HetCount = len(diplos) - (Hom1Count + Hom2Count)
    if verbose: sys.stderr.write(str(Hom1Count) + " " + str(Hom2Count) + " " + str(HetCount))
    p = HWEtest(HetCount,Hom1Count,Hom2Count)
    if verbose: sys.stderr.write("P: " + str(p))
    if p <= P_value: return False
    else: return True


def siteTest(site,samples=None,minCalls=1,minPopCalls=None,minAlleles=0,maxAlleles=float("inf"),
             minPopAlleles=None,maxPopAlleles=None,minVarCount=None,maxHet=None,minFreq=None,maxFreq=None,
             HWE_P=None,HWE_side="both",fixed=False,nearlyFixedDiff=None):
    if not samples: samples = site.sampleNames
    #check sufficient number of non-N calls
    if site.nonMissing() < minCalls: return False
    numBases = site.asList(mode = "numeric", samples=samples)
    numBases = numBases[numBases >= 0]
    #check min and max alleles 
    nAlleles = len(set(site.alleles(samples)))
    if not minAlleles <= nAlleles <= maxAlleles: return False
    #check variant filters
    if nAlleles > 1:
        # minor allele count
        if minVarCount and sorted(binBaseFreqs(numBases, asCounts = True))[-2] < minVarCount: return False
        #check maximum heterozygots?
        if maxHet and site.hets(samples) > maxHet: return False
        #if there is a frequency cutoff
        if minFreq and not minFreq <= sorted(binBaseFreqs(numBases))[-2]: return False
        if maxFreq and not sorted(binBaseFreqs(numBases))[-2] <= maxFreq: return False
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
        for popName in popNames:
            if minPopCalls:
                popCalls = sum([site.genotypes[sample].isMissing()==False for sample in site.pops[popName]])
                if popCalls < minPopCalls[popName]: return False
        #if we want fixed differences only and there are two or more pops specified
        if fixed or minPopAlleles or maxPopAlleles:
            #if fixed all pops must have only one allele, but taken together must have more than one
            allelesByPop = [site.alleles(pop=popName) for popName in popNames]
            if fixed and not (len(set(map(len, allelesByPop))) == 1 and
                              len(set([a for popAlleles in allelesByPop for a in popAlleles]))) > 1: return False
            
            #otherwise 
            if minPopAlleles or maxPopAlleles:
                if minPopAlleles == None: minPopAlleles = dict(zip(popNames, [0]*len(popNames)))
                if maxPopAlleles == None: maxPopAlleles = dict(zip(popNames, [4]*len(popNames)))
                for x in range(len(popNames)):
                    if not minPopAlleles[popNames[x]] <= len(allelesByPop[x]) <= maxPopAlleles[popNames[x]]: return False
        
        #if we want nearly fixed differences, we need to get pop freqs and find any freq difference big enough
        elif nearlyFixedDiff is not None:
            popFreqs = [site.baseFreqs(pop=popName) for popName in popNames]
            freqDiffs = [popFreqs[c[0]] - popFreqs[c[1]] for c in list(itertools.combinations(range(len(popNames)), 2))]
            if not np.any(np.absolute(np.concatenate(freqDiffs)) >= nearlyFixedDiff): return False
    
    #if we get here we've passed all filters
    return True



######################################################################################################################

#modules for working with and analysing alignments


class Alignment:
    __slots__ = ['sequences', 'names', 'groups', 'groupIndDict', 'indGroupDict',
                 'length', 'numArray', 'positions', 'sampleNames', 'nanMask', 'N', 'l',
                 'array','numArray', '_distMat_']
    
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
         
        self.nanMask = self.numArray >= 0
        
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
        
        #we make a None object for distance matrix, but if any dist matrix type function is run, this will be filled
        self._distMat_ = None
    
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
    
    def pairDist(self, i, j):
        nanMask = self.nanMask[i,:] & self.nanMask[j,:]
        return numHamming(self.numArray[i,:][nanMask], self.numArray[j,:][nanMask])
    
    def distMatrix(self):
        distMat = np.zeros((self.N,self.N))
        for i in range(self.N - 1):
            for j in range(i + 1, self.N):
                distMat[i,j] = distMat[j,i] = self.pairDist(i,j)
        self._distMat_ = distMat
        return distMat
    
    def sampleHet(self, sampleNames=None, asList = False):
        if sampleNames is None: sampleNames,sampleIndices = uniqueIndices(self.sampleNames, preserveOrder=True)
        else: sampleIndices = [np.where(self.sampleNames == sampleName)[0] for sampleName in sampleNames]
        #if a pre-computed distance matrix is available, use that
        if self._distMat_ is not None:
            hets = [self._distMat_[x[0],x[1]] if len(x)==2 else np.NaN for x in sampleIndices]
        #otherwise compute pairwise distances (i.e. entire matrix not needed)
        else:
            hets = [self.pairDist(x[0],x[1]) if len(x)==2 else np.NaN for x in sampleIndices]
        return dict(zip(sampleNames,hets)) if not asList else hets
    
    #makes a dict of average distance among samples.
    #if all are haploid, this is just a dictionary of the output of distMatrix()
    #if some have ploidy > 1, this will average distance among sample haplotypes
    def indPairDists(self, asDict=True, includeSameWithSame=False):
        distMat = self.distMatrix() if self._distMat_ is None else self._distMat_ 
        #mask diagonal if necessary
        if not includeSameWithSame: np.fill_diagonal(distMat, np.NaN) # set all same-with-same to Na
        sampleNames,sampleIndices = uniqueIndices(self.sampleNames, preserveOrder=True)
        n = len(sampleNames)
        if asDict:
            pairDists = {}
            for sampleName in sampleNames: pairDists[sampleName] = {} 
            for i,j in itertools.product(range(n),repeat=2):
                pairDists[sampleNames[i]][sampleNames[j]] = np.nanmean(distMat[np.ix_(sampleIndices[i],sampleIndices[j])])
            
            return pairDists
        else:
            indDistMat = np.zeros(n,n)
            for i,j in itertools.combinations_with_replacement(range(n),2):
                    indDistMat[i,j] = indDistMat[j,i] = np.nanmean(distMat[np.ix_(sampleIndices[i],sampleIndices[j])])
            return indDistMat
    
    def groupDistStats(self, doPairs = True):
        #get distance matrix unless a precomputed one is available
        distMat = self._distMat_ if self._distMat_ is not None else self.distMatrix()
        np.fill_diagonal(distMat, np.NaN) # set all same-with-same to Na
        
        pops,indices = np.unique(self.groups, return_inverse = True)
        nPops = len(pops)
        
        #get population indices - which positions in the alignment correspond to each population
        # this will allow indexing specific pops from the matrix.
        popIndices = [list(np.where(indices==x)[0]) for x in range(nPops)]
        
        output = {}
        
        #pi for each pop
        for x in range(nPops):
            output["pi_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[x])])
        
        if nPops == 1 or not doPairs: return output

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
    
    def varSites(self, indices=None, names=None):
        if names is not None: indices = np.where(np.in1d(aln.names,names))[0]
        if indices is None: indices = np.arange(self.N)
        return np.where([numVar(self.numArray[indices,x][self.nanMask[indices,x]]) for x in range(self.l)])[0]
    
    def groupFreqStats(self):
        #dictionary of popgen statistics based on sites (as opposed to pairwise sequence comparisons)
        # THIS ONLY USES SITES WITHOUT ANY MISSING DATA IN A GIVEN GROUP
        output = {}
        
        for groupName in np.unique(self.groups):
            seqIdx = np.where(self.groups==groupName)[0]
            N = len(seqIdx)
            siteIdx = np.where(np.all(self.nanMask, axis=0))[0]
            l = len(siteIdx)
            if l >= 1:
                siteBaseCounts = np.apply_along_axis(binBaseFreqs,0,self.numArray[np.ix_(seqIdx,siteIdx)],asCounts=True)
                #Here I caculate a site-wise pi by summing multplied pairs of base frequencies and dividing by total possible pairs
                #this is slower than computing from the SFS, but it allows for more than 2 alleles per site
                sitePi = np.apply_along_axis(baseCountPi, 0, siteBaseCounts)
                S = sum(sitePi != 0.)
                thetaPi = sum(sitePi)
                thetaW = S / sum(1./np.arange(1,N))
                TajD = TajimaD(N, S, thetaPi)
            else: S=thetaPi=thetaW=TajD=np.NaN
            output["l_"+groupName] = l
            output["S_"+groupName] = S
            output["thetaPi_"+groupName] = thetaPi
            output["thetaW_"+groupName] = thetaW
            output["TajD_"+groupName] = TajD
        
        return output
    
    def biSites(self): return np.where([len(np.unique(self.numArray[:,x][self.nanMask[:,x]])) == 2 for x in range(self.l)])[0]
        
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
    
    def consensus(self, minData = 0.001):
        propData = self.siteNonNan(prop=True)
        return [consensus(self.array[:,i][self.nanMask[:,i]]) if propData[i] >= minData else "N" for i in range(self.l)]
    
    def alleles(self):
        return [set(self.array[:,i][self.nanMask[:,i]]) for i in range(self.l)]
    
    def sampleAlleles(self, sampleNames=None, asList = False):
        if sampleNames is None: sampleNames,sampleIndices = uniqueIndices(self.sampleNames, preserveOrder=True)
        else: sampleIndices = [np.where(self.sampleNames == sampleName)[0] for sampleName in sampleNames]
        #get alleles at each site for each sample
        sampleAlleles = [[set(self.array[sidx,i][self.nanMask[sidx,i]]) for sidx in sampleIndices] for i in range(self.l)]
        if asList: return sampleAlleles
        else: return [dict(zip(sampleNames,sampleAlleles[i])) for i in range(self.l)]


def genoToAlignment(seqDict, sampleData=None, genoFormat = "diplo", positions = None):
    if sampleData is None: sampleData = SampleData()
    seqNames = []
    sampleNames = []
    groups = []
    haploidSeqs = []
    #first pseudo phase all seqs if necessary
    for indName in seqDict.keys():
        seqList = splitSeq(seqDict[indName], genoFormat)
        ploidy = sampleData.ploidy[indName] if indName in sampleData.ploidy and sampleData.ploidy[indName] != None else len(seqList)
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
                     positions= positions,
                     sampleNames= [sampleNames[i] for i in order])




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
    #for x in range(len(seqA)):
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
    #get distance matrix unless a precomputed one is available
    distMat = Aln._distMat_ if Aln._distMat_ is not None else Aln.distMatrix()
    np.fill_diagonal(distMat, np.NaN) # set all same-with-same to Na
    
    pops,indices = np.unique(Aln.groups, return_inverse = True)
    nPops = len(pops)
    
    #get population indices - which positions in the alignment correspond to each population
    # this will allow indexing specific pops from the matrix.
    popIndices = [list(np.where(indices==x)[0]) for x in range(nPops)]
    
    output = {}
    
    #pi for each pop
    for x in range(nPops):
        output["pi_" + pops[x]] = np.nanmean(distMat[np.ix_(popIndices[x],popIndices[x])])
    
    if nPops == 1 or not doPairs: return output

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
    #ABBA = BABA = D_numer = D_denom = fd_denom = fdM_denom = 0.0
    #sitesUsed = 0
    ##get derived frequencies for all biallelic siites
    #for i in P123Aln.biSites():
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
        ## get weigtings for ABBAs and BABAs
        #_ABBA_ = (1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
        #_BABA_ = P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
        #ABBA += _ABBA_
        #BABA += _BABA_
        #D_numer += _ABBA_ - _BABA_
        #D_denom += _ABBA_ + _BABA_
        
        ##fd
        #if P2derFreq <= P3derFreq:
            #fd_denom += (1 - P1derFreq) * P3derFreq * P3derFreq * (1 - P4derFreq) - P1derFreq * (1 - P3derFreq) * P3derFreq * (1 - P4derFreq)
        #else:
            #fd_denom += (1 - P1derFreq) * P2derFreq * P2derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P2derFreq * (1 - P4derFreq)
        
        ##fdM
        #if P2derFreq >= P1derFreq:
            #if P2derFreq <= P3derFreq:
                #fdM_denom += (1 - P1derFreq) * P3derFreq * P3derFreq * (1 - P4derFreq) - P1derFreq * (1 - P3derFreq) * P3derFreq * (1 - P4derFreq)
            #else:
                #fdM_denom += (1 - P1derFreq) * P2derFreq * P2derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P2derFreq * (1 - P4derFreq)
        #else:
            #if P1derFreq <= P3derFreq:
                #fdM_denom -= (1 - P3derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq) - P3derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
            #else:
                #fdM_denom -= (1 - P1derFreq) * P2derFreq * P1derFreq * (1 - P4derFreq) - P1derFreq * (1 - P2derFreq) * P1derFreq * (1 - P4derFreq)
        
        #sitesUsed += 1
    ##calculate D, fd
    #output = {}
    #output["ABBA"] = ABBA
    #output["BABA"] = BABA
    #try: output["D"] = D_numer / D_denom
    #except: output["D"] = np.NaN
    #try:
        #if output["D"] >= 0: output["fd"] = D_numer / fd_denom
        #else: output["fd"] = np.NaN
    #except: output["fd"] = np.NaN
    #try:
        #output["fdM"] = D_numer / fdM_denom
    #except: output["fd"] = np.NaN
    #output["sitesUsed"] = sitesUsed
    
    #return output

### F4 and ABBABABA stats using Hannes Svardal's method

def f4(p1,p2,p3,p4):
    f4 = (1 - p1)*p2*p3*(1-p4) - p1 * (1-p2)*p3*(1-p4)
    return f4

def f4_c(p1,p2,p3,p4):
    """
    corrected f4
    """
    f = f4(p1,p2,p3,p4) + f4(1-p1,1-p2,1-p3,1-p4)
    return f

def fhom_old(p1,p2,p3,p4):
    return f4(p1,p2,p3,p4).sum()*1./f4(p1,p3,p3,p4).sum()

def fhom_new(p1,p2,p3,p4):
    """
    This is fhom using the correct formula
    for f4
    """
    return (f4_c(p1,p2,p3,p4)).sum()*1./(f4_c(p1,p3,p3,p4)).sum()

def D(p1,p2,p3,p4):
    return f4(p1,p2,p3,p4).sum()*1./((1 - p1)*p2*p3*(1-p4) + p1 * (1-p2)*p3*(1-p4)).sum()



def D_new(p1,p2,p3,p4):
    """
    This is D using the correct formula
    for f4.
    One benchmark for a correct statistic is that
    it should always have the same sign as D_new!
    """
    return (f4_c(p1,p2,p3,p4)).sum()*1./((1 - p1)*p2*p3*(1-p4) + p1 * (1-p2)*p3*(1-p4)+\
                                                              p1*(1-p2)*(1-p3)*p4 + (1-p1) * p2*(1-p3)*p4).sum()

def fd(p1,p2,p3,p4):
    pd = p2* (p2>p3) + p3*(p3>=p2)
    return f4(p1,p2,p3,p4).sum()*1./f4(p1,pd,pd,p4).sum()


def fd_new(p1,p2,p3,p4):
    """
    This is fd using the correct formula
    for f4.
    I think it generally does not make sense 
    with the general f4.
    """
    pd = p2* (p2>p3) + p3*(p3>=p2)
    return (f4_c(p1,p2,p3,p4)).sum()*1./(f4_c(p1,pd,pd,p4)).sum()

def get_fdm_p(p1,p2,p3):
    a = (p3 > p1)
    b = (p3 > p2)
    x = (p1 > p2)
    y = ~x
    pdm1 = p3*(x&a) + p1*(~(x&a))
    pdm2 = p3*(y&b) + p2*(~(y&b))
    pdm3 = -p3*(x&a) + p3*(y&b) - p1*(x&~a) + p2*(y&~b)
    return pdm1, pdm2, pdm3

def fdm(p1,p2,p3,p4):
    pdm1,pdm2,pdm3 = get_fdm_p(p1,p2,p3)
    denom = f4(pdm1,pdm2,pdm3,p4)
    #This would be a new version wiht the correct f formula, but I think it does not make much sense
    #denom = ((pdm1 - pdm2) * (pdm3 - p4))
    return f4(p1,p2,p3,p4).sum()*1./denom.sum()


def fdm_new(p1,p2,p3,p4):
    """
    This is fd using the correct formula
    for f4.
    I think it generally does not make sense 
    with the general f4.
    """
    pdm1,pdm2,pdm3 = get_fdm_p(p1,p2,p3)
    #denom = ((pdm1 - pdm2) * (pdm3 - p4))
    denom = f4_c(pdm1,pdm2,pdm3,p4)
    return (f4_c(p1,p2,p3,p4)).sum()*1./denom.sum()

def fdh(p1,p2,p3,p4):
    """
    This is fd using the correct formula
    for f4.
    I think it generally does not make sense 
    with the general f4.
    """
    num = f4_c(p1,p2,p3,p4)
    t11 = f4_c(p1,p3,p3,p4) 
    t12 = f4_c(p4,p2,p3,p4)
    t21 = f4_c(p3,p2,p3,p4)
    t22 = f4_c(p1,p4,p3,p4)
    denom = np.amax([t11,t12,t21,t22],axis=0)
    
    return num.sum()*1./denom.sum()

def fdh2(p1,p2,p3,p4):
    """
    This is fd using the correct formula
    for f4.
    I think it generally does not make sense 
    with the general f4.
    """
    num = f4_c(p1,p2,p3,p4)
    t11 = f4_c(p1,p3,p3,p4) 
    t12 = f4_c(p4,p2,p3,p4)
    t21 = f4_c(p3,p2,p3,p4)
    t22 = f4_c(p1,p4,p3,p4)
    t31 = f4_c(p1,p2,p2,p4)
    t32 = f4_c(p1,p2,p3,p1)
    t41 = f4_c(p1,p2,p1,p4)
    t42 = f4_c(p1,p2,p3,p2)
    denom = np.amax([t11,t12,t21,t22,t31,t32,t41,t42],axis=0)
    
    return num.sum()*1./denom.sum()



def fh(p1,p2,p3,p4):
    """
    PROPOSED statistic
    I think that this statistic is well behaved
    in the sense that it is in [-1,1] and
    it always has the same sign as D_new.
    I am not sure whether it is nicely proportional to the amount of gene flow.
    Actually, I think that one cannot have both an unbiased estimator of gene flow
    and totally nice behaviour even for small numbers.
    I can imagine any esimator that reduced variability (with respect to f_green)
    automatically introduces a bias.
    
    This proposed f is basically (assuming the correct f4):
    
    sum(f4(p1,p2,p3,p4))
    --------------------
    sum(max(f4(p3,p4,p3,p4),f4(p1,p2,p1,p2)))
    
    or, equivalently
    
    sum(f4(p1,p2,p3,p4))
    --------------------
    sum(max(f2(p3,p4),f2(p1,p2)))
    
    """
    t1 = np.abs((p1-p2))
    t2 = np.abs((p3-p4))
    denom = (t1 * (t1>t2) + t2 * (t2>=t1))**2
    return (f4(p1,p2,p3,p4)+f4(1-p1,1-p2,1-p3,1-p4)).sum()*1./denom.sum()


def ABAA(p1,p2,p3,p4):
    return ((1 - p1)*p2*(1-p3)*(1-p4)).sum()

def BAAA(p1,p2,p3,p4):
    return (p1*(1 - p2)*(1-p3)*(1-p4)).sum()

def ABBA(p1,p2,p3,p4):
    return ((1 - p1)*p2*p3*(1-p4)).sum()

def BABA(p1,p2,p3,p4):
    return (p1*(1-p2)*p3*(1-p4)).sum()

def ABAA_BABB(p1,p2,p3,p4):
    return ((1 - p1)*p2*(1-p3)*(1-p4) + p1*(1-p2)*p3*p4).sum()

def BAAA_ABBB(p1,p2,p3,p4):
    return (p1*(1 - p2)*(1-p3)*(1-p4) + (1-p1)*p2*p3*p4).sum()

def ABBA_BAAB(p1,p2,p3,p4):
    return ((1 - p1)*p2*p3*(1-p4) + p1*(1-p2)*(1-p3)*p4).sum()

def BABA_ABAB(p1,p2,p3,p4):
    return (p1*(1-p2)*p3*(1-p4) + (1-p1)*p2*(1-p3)*p4).sum()


##new fourPop (previously ABBABABA) code, implementing all of Hannes Svardal's methods, plus others based on numpy arrays
def fourPop(aln, P1, P2, P3, P4, minData, polarize=False, fixed=False):
    #subset by population
    all4Aln = aln.subset(groups=[P1,P2,P3,P4])
    P1Aln = all4Aln.subset(groups=[P1])
    P2Aln = all4Aln.subset(groups=[P2])
    P3Aln = all4Aln.subset(groups=[P3])
    P4Aln = all4Aln.subset(groups=[P4])

    biallelic = [len(np.unique(all4Aln.numArray[:,x][all4Aln.nanMask[:,x]])) == 2 for x in range(all4Aln.l)]
    
    enoughData =((P1Aln.siteNonNan()*1./P1Aln.N >= minData) &
                    (P2Aln.siteNonNan()*1./P2Aln.N >= minData) &
                    (P3Aln.siteNonNan()*1./P3Aln.N >= minData) &
                    (P4Aln.siteNonNan()*1./P4Aln.N >= minData))
        
    goodSites = np.where(biallelic & enoughData)[0]

    all4freqs = all4Aln.siteFreqs(sites=goodSites)
    P1freqs = P1Aln.siteFreqs(sites=goodSites)
    P2freqs = P2Aln.siteFreqs(sites=goodSites)
    P3freqs = P3Aln.siteFreqs(sites=goodSites)
    P4freqs = P4Aln.siteFreqs(sites=goodSites)
    
    #try:
    if len(goodSites) >= 1:
        if polarize: alleleIndex = np.where((all4freqs > 0) & (P4freqs == 0))
        elif fixed: alleleIndex = np.where((all4freqs > 0) & (P4freqs == 0) &
                                        ((P1freqs==0) | (P1freqs==1)) &
                                        ((P2freqs==0) | (P2freqs==1)) & 
                                        ((P3freqs==0) | (P3freqs==1)))
        else: alleleIndex = (np.arange(all4freqs.shape[0]), np.argsort(all4freqs, axis = 1)[:,2],)
        
        p1 = P1freqs[alleleIndex[0],alleleIndex[1]]
        p2 = P2freqs[alleleIndex[0],alleleIndex[1]]
        p3 = P3freqs[alleleIndex[0],alleleIndex[1]]
        p4 = P4freqs[alleleIndex[0],alleleIndex[1]]
        
        f1 = fhom_old(p1,p2,p3,p4)
        f2 = fhom_new(p1,p2,p3,p4)
        d = D(p1,p2,p3,p4)
        d_new = D_new(p1,p2,p3,p4)
        fd1 = fd(p1,p2,p3,p4)
        fd_new1 = fd_new(p1,p2,p3,p4)
        fdm1 = fdm(p1,p2,p3,p4)
        fdm_new1 = fdm_new(p1,p2,p3,p4)
        fdh1 = fdh(p1,p2,p3,p4)
        fdh21 = fdh2(p1,p2,p3,p4)
        fh1 = fh(p1,p2,p3,p4)
        abba = ABBA(p1,p2,p3,p4)
        baba = BABA(p1,p2,p3,p4)
        abaa = ABAA(p1,p2,p3,p4)
        baaa = BAAA(p1,p2,p3,p4)
        sitesUsed = len(alleleIndex[0])
        
        return dict(zip(['fhom',"fhom'",'D','fd',"fd'",'fdm',"fdm'",'fdh','fdh2','fh',"ABBA","BABA","ABAA","BAAA","sitesUsed"],
                        [f1,f2,d,fd1,fd_new1,fdm1,fdm_new1,fdh1,fdh21,fh1,abba,baba,abaa,baaa,sitesUsed]))
    else:
        return dict(zip(['fhom',"fhom'",'D','fd',"fd'",'fdm',"fdm'",'fdh','fdh2','fh',"ABBA","BABA","ABAA","BAAA","sitesUsed"],
                        [np.NaN]*14 +[0]))


##new ABBABABA code using numpy arrays
def ABBABABA(aln, P1, P2, P3, P4, minData, polarize=True, fixed=False):
    #subset by population
    all4Aln = aln.subset(groups=[P1,P2,P3,P4])
    P1Aln = all4Aln.subset(groups=[P1])
    P2Aln = all4Aln.subset(groups=[P2])
    P3Aln = all4Aln.subset(groups=[P3])
    P4Aln = all4Aln.subset(groups=[P4])

    biallelic = [len(np.unique(all4Aln.numArray[:,x][all4Aln.nanMask[:,x]])) == 2 for x in range(all4Aln.l)]
    
    enoughData =((P1Aln.siteNonNan()*1./P1Aln.N >= minData) &
                    (P2Aln.siteNonNan()*1./P2Aln.N >= minData) &
                    (P3Aln.siteNonNan()*1./P3Aln.N >= minData) &
                    (P4Aln.siteNonNan()*1./P4Aln.N >= minData))
        
    goodSites = np.where(biallelic & enoughData)[0]

    all4freqs = all4Aln.siteFreqs(sites=goodSites)
    P1freqs = P1Aln.siteFreqs(sites=goodSites)
    P2freqs = P2Aln.siteFreqs(sites=goodSites)
    P3freqs = P3Aln.siteFreqs(sites=goodSites)
    P4freqs = P4Aln.siteFreqs(sites=goodSites)
    
    #try:
    if len(goodSites) >= 1:
        if polarize: alleleIndex = np.where((all4freqs > 0) & (P4freqs == 0))
        elif fixed: alleleIndex = np.where((all4freqs > 0) & (P4freqs == 0) &
                                        ((P1freqs==0) | (P1freqs==1)) &
                                        ((P2freqs==0) | (P2freqs==1)) & 
                                        ((P3freqs==0) | (P3freqs==1)))
        else: alleleIndex = (np.arange(all4freqs.shape[0]), np.argsort(all4freqs, axis = 1)[:,2],)
        
        p1 = P1freqs[alleleIndex[0],alleleIndex[1]]
        p2 = P2freqs[alleleIndex[0],alleleIndex[1]]
        p3 = P3freqs[alleleIndex[0],alleleIndex[1]]
        p4 = P4freqs[alleleIndex[0],alleleIndex[1]]
        
        _D_ = D(p1,p2,p3,p4)
        _fd_ = fd(p1,p2,p3,p4)
        _fdm_ = fdm(p1,p2,p3,p4)
        _ABBA_ = ABBA(p1,p2,p3,p4)
        _BABA_ = BABA(p1,p2,p3,p4)
        sitesUsed = len(alleleIndex[0])
        
        return dict(zip(['D','fd','fdM',"ABBA","BABA","sitesUsed"],
                        [_D_,_fd_,_fdm_,_ABBA_,_BABA_,sitesUsed]))
    else:
        return dict(zip(['D','fd','fdM',"ABBA","BABA","sitesUsed"],
                        [np.NaN]*6 +[0]))


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
    __slots__ = ['scaffold', 'limits','sites', 'names', 'positions', 'ID', 'n']
    
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
        end = len(self.positions)
        while i < end and self.positions[i] < self.limits[0]: i += 1
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



#new test implementation of GenoWindow using deque 23 Oct 2018

#from collections import deque

#class GenoWindow:
    #def __init__(self, scaffold = None, limits=[-np.inf,np.inf],sites = None, names = None, positions = None, ID = None):
        #self.names = names if names is not None else []
        #self.n = len(self.names)
        #self.scaffold = scaffold
        #self.limits = limits
        #self.ID = ID
        #if sites is not None and positions is not None:
            #if not len(sites) == len(positions) == 0:
                #assert len(set([len(site) for site in sites])) == 1, "Number of genotypes per site must be equal."
                #assert len(sites[0]) == len(names), "Number of names must match number of genotypes per site."
                #assert len(positions) == len(sites), "Positions must match number of sites"
        #else:
            #sites = deque(maxlen=self.limits[1]-self.limits[0]+1 if self.limits[0] > -np.inf and self.limits[1] < np.inf else None)
            #positions = deque(maxlen=self.limits[1]-self.limits[0]+1 if self.limits[0] > -np.inf and self.limits[1] < np.inf else None)
        #self.sites = sites if sites is not None else deque()
        #self.positions = positions if positions is not None else deque()
    
    #def copy(self): return GenoWindow(scaffold=self.scaffold, limits=self.limits[:],
                                      #sites=copy(self.sites), names=self.names[:], positions=copy(self.positions), ID=self.ID)
    
    ###method for adding
    #def addBlock(self, sites, positions):
        #assert len(set([len(site) for site in sites])) == 1, "Number of genotypes per site must be equal."
        #assert len(sites[0]) == self.n, "Number of genotypes per site must match number of names."
        #assert len(positions) == len(sites), "Positions must match number of sites"
        #assert np.all(self.limits[0] <= positions <= self.limits[1]), "Position outside of window limit"
        #self.sites += sites
        #self.positions += positions
    
    #def addSite(self, GTs, position=np.NaN, ignorePosition=False):
        #assert len(GTs) == self.n, "Number of genotypes per site must match number of names."
        #if not ignorePosition:
            #assert self.limits[0] <= position <= self.limits[1], "Position: " + str(position) + " outside of window limits: " + "-".join([str(l) for l in self.limits])
        #else: position = np.NaN
        #self.positions.append(position)
        #self.sites.append(GTs)
    
    #def seqLen(self): return len(self.positions)
    
    #def firstPos(self): return min(self.positions)
    
    #def lastPos(self): return max(self.positions)
    
    #def slide(self,step=None,newLimits=None):
        ##function to slide window along scaffold
        #assert (step != None and step >= 1) or (newLimits != None and newLimits[0] > self.limits[0]), "Step must be > 0" 
        #if step: self.limits = [l+step for l in self.limits]
        #else: self.limits = newLimits
        #if len(self.positions) == 0: return
        #i = self.positions.popleft()
        #self.sites.popleft()
        #while i < self.limits[0]:
            #try:
                #i = self.positions.popleft()
                #self.sites.popleft()
            #except:
                #break
    
    #def trim(self,right=False,remove=None,leave=None):
        #assert remove != None or leave != None
        #if not remove: remove=self.seqLen() - leave
        #if not right:
            ##trim positions and sites
            #self.positions = self.positions[remove:]
            #self.sites = self.sites[remove:]
        #else:
            #self.positions = self.positions[:-remove]
            #self.sites = self.sites[:-remove]
    
    #def seqDict(self, names=None):
        #if names is None: names = self.names
        #indices = [self.names.index(n) for n in names]
        #return dict(zip(names, [[site[i] for site in self.sites] for i in indices]))
    
    #def midPos(self):
        #try: return int(round(sum(self.positions)/len(self.positions)))
        #except: return np.NaN

def parseGenoLine(line, names, scafCol=0, posCol=1, firstSampleCol=2,
                  type=str, splitPhased=False, asDict = True, precompDict=None, addToPrecomp=True):
    if line and line != "":
        lineData = line.split(None,firstSampleCol)
        GTstring = lineData[-1]
        #check if there is a precompiled genotype dictionary, and if so, check for the line in there
        if precompDict and GTstring in precompDict:
            GTs = precompDict[GTstring]
        else:
            GTs = GTstring.split()
            if splitPhased: GTs = [a for GT in GTs for a in list(GT)[::2]]
            if type!=str: GTs = [float(GT) if type==float else int(GT) for GT in GTs]
            if asDict: GTs = dict(zip(names,GTs))
            if (precompDict != None) and addToPrecomp:
                precompDict[GTstring] = GTs
                precompDict["__counter__"] += 1
        return {"scaffold": lineData[scafCol] if scafCol >= 0 else None,
                "position": int(lineData[posCol]) if posCol >= 0 else None,
                "GTs": GTs}
    else:
        return {"scaffold": None, "position": None, "GTs": None}



#a little function to get the next from a generator while first checking pythoin version
#will be redundant when completely ported to python 3

def getNext(generator):
    return generator.next() if sys.version_info.major < 3 else next(generator)

class GenoFileReader:
    def __init__(self,genoFile, headerLine=None, scafCol=0, posCol=1, firstSampleCol=2,
                 type=str, splitPhased = False, ploidy=None, precomp=True, precompMaxSize=10000):
        self.genoFile = genoFile
        if not headerLine: headerLine = getNext(genoFile)
        self.names = headerLine.split()[firstSampleCol:]
        self.scafCol = scafCol
        self.posCol = posCol
        self.firstSampleCol = firstSampleCol
        self.type = type
        self.splitPhased=splitPhased
        if splitPhased:
            assert ploidy is not None, "Ploidy must be defined for splitting phased sequences"
            if self.names:
                self.names = makeHaploidNames(self.names, ploidy)
        #add a dictionary for precompiled genotypes, if you want one
        self.precompDict = {}
        self.precompDict["__maxSize__"] = precompMaxSize
        self.precompDict["__counter__"] = 0
    
    def siteBySite(self, asDict=True):
        for line in self.genoFile:
            yield parseGenoLine(line, self.names, self.scafCol, self.posCol, self.firstSampleCol, self.type, self.splitPhased, asDict,
                                self.precompDict, addToPrecomp=self.precompDict["__counter__"]<self.precompDict["__maxSize__"])
    
    def nextSite(self, asDict=True):
        try: line = getNext(self.genoFile)
        except: line = None
        return parseGenoLine(line, self.names, self.scafCol, self.posCol, self.firstSampleCol, self.type, self.splitPhased, asDict,
                            self.precompDict, addToPrecomp=self.precompDict["__counter__"]<self.precompDict["__maxSize__"])


#function to read entire genoFile into a window-like object
def parseGenoFile(genoFile, headerLine=None, names = None, includePositions = False, splitPhased=False, ploidy=None):
    #file reader
    reader=GenoFileReader(genoFile, headerLine, splitPhased=splitPhased, ploidy=ploidy)
    #get names (only needed if we don't want to read all sequences in the file, otherwise we just get them from the file)
    if names:
        extractSpecificGTs = True
        if splitPhased: names = makeHaploidNames(names, ploidy)
    else:
        extractSpecificGTs = False
        names = reader.names
    
    #initialise window
    window = GenoWindow(names = names)
    #populate window
    for siteData in reader.siteBySite(asDict=extractSpecificGTs):
        GTs = [siteData["GTs"][name] for name in names] if extractSpecificGTs else siteData["GTs"]
        window.addSite(GTs=GTs, position=siteData["position"], ignorePosition= not includePositions)
    
    return window


#sliding window generator function
def slidingCoordWindows(genoFile, windSize, stepSize, names = None, splitPhased=False, ploidy = None,
                        include = None, exclude = None, skipDeepcopy = False):
    #file reader
    reader=GenoFileReader(genoFile, splitPhased=splitPhased, ploidy=ploidy)
    #get names
    if names:
        extractSpecificGTs = True
        if splitPhased: names = makeHaploidNames(names, ploidy)
    else:
        extractSpecificGTs = False
        names = reader.names
    #window counter
    windowsDone = 0
    #initialise window
    window = GenoWindow(names = names)
    #first site
    site = reader.nextSite(asDict = extractSpecificGTs)
    while site["position"] is not None:
        #build window
        while site["scaffold"] == window.scaffold and site["position"] <= window.limits[1]:
            if site["position"] >= window.limits[0]:
                #add this site to the window
                GTs = [site["GTs"][name] for name in names] if extractSpecificGTs else site["GTs"]
                window.addSite(GTs=GTs, position=site["position"])
            #read next line
            site = reader.nextSite(asDict = extractSpecificGTs)
        
        '''if we get here, the line in hand is incompatible with the currrent window
            If the window is not empty, yield it'''
        
        if window.scaffold is not None:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield window.copy()
        
        #now we need to make a new window
        #if on same scaffold, just slide along
        if site["scaffold"] == window.scaffold:
            window.slide(step = stepSize)
            window.ID = windowsDone + 1
        
        #otherwise we're on a new scaffold (or its the end of the file)
        else:
            #if its one we want to analyse, start new window
            if (not include and not exclude) or (include and site["scaffold"] in include) or (exclude and site["scaffold"] not in exclude):
                window = GenoWindow(scaffold = site["scaffold"], limits=[1,windSize], names = names, ID = windowsDone + 1)
            
            #if its a scaf we don't want, were going to read lines until we're on one we do want
            else:
                badScaf = site["scaffold"]
                while site["scaffold"] == badScaf or (include and site["scaffold"] not in include and site["scaffold"] is not None) or (exclude and site["scaffold"] in exclude and site["scaffold"] is not None):
                    site = reader.nextSite(asDict = extractSpecificGTs)
            
        #if we've reached the end of the file, break
        if site["position"] is None:
            break
    


#sliding window generator function
def slidingSitesWindows(genoFile, windSites, overlap, maxDist = np.inf, minSites = None, names = None,
                        splitPhased=False, ploidy=None, include = None, exclude = None, skipDeepcopy = False):
    if not minSites: minSites = windSites #if minSites < eindSites, windows at ends of scaffolds can still be emmitted
    #file reader
    reader=GenoFileReader(genoFile, splitPhased=splitPhased, ploidy=ploidy)
    #get names
    if names:
        extractSpecificGTs = True
        if splitPhased: names = makeHaploidNames(names, ploidy)
    else:
        extractSpecificGTs = False
        names = reader.names
    #window counter
    windowsDone = 0
    #initialise window
    window = GenoWindow(names = names)
    #first site
    site = reader.nextSite(asDict = extractSpecificGTs)
    while site:
        #build window
        while site["scaffold"] == window.scaffold and window.seqLen() < windSites and (window.seqLen() == 0 or site["position"] - window.firstPos() <= maxDist):
            #add this site to the window
            GTs = [site["GTs"][name] for name in names] if extractSpecificGTs else site["GTs"]
            window.addSite(GTs=GTs, position=site["position"])
            #read next line
            site = reader.nextSite(asDict = extractSpecificGTs)
        
        '''if we get here, either the window is full, or the line in hand is incompatible with the currrent window
            If the window has more than minSites, yield it'''
        
        
        if window.seqLen() >= minSites:
            windowsDone += 1
            
            if skipDeepcopy: yield window
            else: yield deepcopy(window)
            
            #now we need to make a new window
            #if on same scaffold, just trim
            if site["scaffold"] == window.scaffold:
                window.trim(leave = overlap)
                window.ID = windowsDone + 1
            
            #otherwise we're on a new scaffold (or its the end of the file)
            else:
                #if its one we want to analyse, start new window
                if (not include and not exclude) or (include and site["scaffold"] in include) or (exclude and site["scaffold"] not in exclude):
                    window = GenoWindow(scaffold = site["scaffold"], names = names, ID = windowsDone + 1)
                
                #if its a scaf we don't want, were going to read lines until we're on one we do want
                else:
                    badScaf = site["scaffold"]
                    while site["scaffold"] == badScaf or (include and site["scaffold"] not in include and site["scaffold"] is not None) or (exclude and site["scaffold"] in exclude and site["scaffold"] is not None):
                    
                        site = reader.nextSite(asDict = extractSpecificGTs)
        
        #If there are insufficient sites, and we're on the same scaffold, just trim off the furthest left site
        else:
            if site["scaffold"] == window.scaffold:
                window.trim(remove = 1)
            
            #If we're on a new scaffold, we do as above
            else:
                #if its one we want to analyse, start new window
                if (not include and not exclude) or (include and site["scaffold"] in include) or (exclude and site["scaffold"] not in exclude):
                    window = GenoWindow(scaffold = site["scaffold"], names = names, ID = windowsDone + 1)
                
                #if its a scaf we don't want, were going to read lines until we're on one we do want
                else:
                    badScaf = site["scaffold"]
                    while site["scaffold"] == badScaf or (include and site["scaffold"] not in include and site["scaffold"] is not None) or (exclude and site["scaffold"] in exclude and site["scaffold"] is not None):
                    
                        site = reader.nextSite(asDict = extractSpecificGTs)
        
        #if we've reached the end of the file, break
        if site["position"] is None:
            break


#window generator function using pre-defined coordinates
def predefinedCoordWindows(genoFile, windCoords, names = None, splitPhased=False, ploidy=None, skipDeepcopy = False):
    #get the order of scaffolds
    allScafs = [w[0] for w in windCoords]
    scafs = sorted(set(allScafs), key=lambda x: allScafs.index(x))
    #file reader
    reader=GenoFileReader(genoFile, splitPhased=splitPhased, ploidy=ploidy)
    #get names
    if names:
        extractSpecificGTs = True
        if splitPhased: names = makeHaploidNames(names, ploidy)
    else:
        extractSpecificGTs = False
        names = reader.names
    window = None
    #read first line
    site = reader.nextSite(asDict = extractSpecificGTs)
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
        while site["scaffold"] and (site["scaffold"] not in scafs or scafs.index(site["scaffold"]) < windScafIdx):
            badScaf = site["scaffold"]
            while site["scaffold"] == badScaf:
                    site = reader.nextSite(asDict = extractSpecificGTs)
        
        #if we're on the right scaffold but abve thwe windiow, keep reading
        while site["scaffold"] == window.scaffold and site["position"] < window.limits[0]:
            site = reader.nextSite(asDict = extractSpecificGTs)
        
        #if we are in a window - build it
        while site["scaffold"] == window.scaffold and window.limits[0] <= site["position"] <= window.limits[1]:
            #add this site to the window
            GTs = [site["GTs"][name] for name in names] if extractSpecificGTs else site["GTs"]
            window.addSite(GTs=GTs, position=site["position"])
            #read next line
            site = reader.nextSite(asDict = extractSpecificGTs)
        
        '''When we get here, either:
            We're on the right scaffold but below the current window
            We're on a scaffold in the list but below the current window
            We've reached the end of the file
            So we have to yield the current window'''
        
        if skipDeepcopy: yield window
        else: yield deepcopy(window)
        
        if site["position"] is None:
            break


#function to read blocks of n lines
#sliding window generator function
def nonOverlappingSitesWindows(genoFile, windSites, names = None, splitPhased=False, ploidy=None, include = None, exclude = None):
    #file reader
    reader=GenoFileReader(genoFile, splitPhased=splitPhased, ploidy=ploidy)
    #get names
    if names:
        extractSpecificGTs = True
        if splitPhased: names = makeHaploidNames(names, ploidy)
    else:
        extractSpecificGTs = False
        names = reader.names
    #blocks counter
    windowsDone = 0
    #initialise an empty block
    window = None
    #read first line
    site = reader.nextSite(asDict = extractSpecificGTs)
    while True:
        #initialise window
        #if its a scaffold we want to analyse, start new window
        if (not include and not exclude) or (include and site["scaffold"] in include) or (exclude and site["scaffold"] not in exclude):
            window = GenoWindow(scaffold = site["scaffold"], names = names, ID = windowsDone + 1)
            
        #if its a scaf we don't want, were going to read lines until we're on one we do want
        else:
            window = None
            badScaf = site["scaffold"]
            while site["scaffold"] == badScaf or (include and site["scaffold"] not in include and site["scaffold"] is not None) or (exclude and site["scaffold"] in exclude and site["scaffold"] is not None):
            
                site = reader.nextSite(asDict = extractSpecificGTs)

        #build window
        while window and site["scaffold"] == window.scaffold and window.seqLen() < windSites:
            #add this site to the window
            GTs = [site["GTs"][name] for name in names] if extractSpecificGTs else site["GTs"]
            window.addSite(GTs=GTs, position=site["position"])
            #read next line
            site = reader.nextSite(asDict = extractSpecificGTs)
        
        '''if we get here, either the window is full, or the line in hand is incompatible with the currrent window
            If the window has more than minSites, yield it'''
        
        if window:
            windowsDone += 1
            yield window
                
        #if we've reached the end of the file, break
        if site["position"] is None:
            break



##########################################################################################################

#functions to make and parse alignment strings in fasta or phylip format


def makeAlnString(names=None, seqs=None, seqDict=None, outFormat="phylip", lineLen=None, NtoGap=False):
    assert outFormat=="phylip" or outFormat=="fasta"
    if seqDict: names, seqs = list(zip(*seqDict.items()))
    else: assert len(names) == len(seqs)
    seqs = ["".join(s) for s in seqs]
    if NtoGap: seqs = [seq.translate(missingtrans) for seq in seqs]
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

#Phylip can include multiple alignments. This will output multiple sequences as a list of tuples
#Single alignments will be output as a tuple, unless asList is True
def parsePhylip(string, asList=False): 
    lineParts = [l.strip().split() for l in string.strip().split("\n")]
    lineParts = [parts for parts in lineParts if parts != []] 
    headIdx = []
    Ns = []
    Ls = []
    for x in range(len(lineParts)):
        try:
            Ls.append(int(lineParts[x][1]))
            Ns.append(int(lineParts[x][0]))
            headIdx.append(x)
        except: pass
    
    headIdx.append(len(lineParts))
    names = [[lineParts[headIdx[i]+1+j][0] for j in range(Ns[i])] for i in range(len(headIdx)-1)]
    seqIdx = [[range(headIdx[i]+1+j,headIdx[i+1],Ns[i]) for j in range(Ns[i])] for i in range(len(headIdx)-1)]
    seqs = [["".join([lineParts[y][1] for y in x]) for x in w] for w in seqIdx] 
    if not asList and len(names) == 1: return (names[0], seqs[0])
    else: return list(zip(names,seqs))


############### writing distance matrices

def makeDistMatString(distArray, roundTo=10):
    return "\n".join([" ".join(i) for i in distArray.round(roundTo).astype(str)])

def makeDistMatPhylipString(distArray, names, roundTo=10):
    output = str(distArray.shape[0]) + "\n"
    for i in range(len(names)):
        output += str(names[i]) + "  " + " ".join(distArray[i,:].round(roundTo).astype(str)) + "\n"
    return output

def makeDistMatNexusString(distArray, names, roundTo=10):
    output = "\nBEGIN Taxa;\nDIMENSIONS ntax={};\nTAXLABELS\n".format(len(names))
    for i in range(len(names)):
        output += "[{}] '{}'\n".format(i+1,names[i])
    output += ";\nEND; [Taxa]\n"
    output += "\nBEGIN Distances;\nDIMENSIONS ntax={};\nFORMAT labels=left diagonal triangle=both;\nMATRIX\n".format(len(names))
    for i in range(len(names)):
        output += "[{}] '{}'    ".format(i+1,names[i]) + " ".join(distArray[i,:].round(roundTo).astype(str)) + "\n"
    output += ";\nEND; [Distances]\n"
    return output


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

###########################################################################################################
#class for working with intervals

def parseRegionText(regionText):
    splitText = regionText.split(":")
    seqName = splitText[0]
    if len(splitText) < 3 or splitText[2] == "": ori = "+"
    else: ori = splitText[2]
    if ori not in "+-": raise ValueError("Incorrect region specification")
    try:
        fromTo = [int(x) for x in splitText[1].split("-")]
        if len(fromTo) == 1: fromTo.append(None)
        if fromTo[1] != None and fromTo[0] > fromTo[1]:
            fromTo = fromTo[::-1]
            ori = "-"
        return (seqName,fromTo[0],fromTo[1],ori,)
    except: return (seqName,None,None,ori,)


def parseRegionList(regionList):
    seqName = regionList[0]
    if len(regionList) < 4: ori = "+"
    else: ori = regionList[3]
    if ori not in "+-": raise ValueError("Orientation must be + or -")
    try:
        fromTo = [int(x) for x in regionList[1:3]]
        if len(fromTo) == 1: fromTo.append(None)
        if fromTo[1] != None and fromTo[0] > fromTo[1]:
            fromTo = fromTo[::-1]
            ori = "-"
        return (seqName,fromTo[0],fromTo[1],ori,)
    except: return (seqName,None,None,ori,)


class Intervals():
    def __init__(self, regions=None, tuples=None, chroms=None, starts=None, ends=None):
        if regions is not None:
            tuples = [parseRegionText(r) for r in regions]
        if tuples is not None:
            self.chroms = np.array([t[0] for t in tuples])
            self.starts = np.array([t[1] if len(t)>1 and t[1] is not None else 0 for t in tuples], dtype=float)
            self.ends = np.array([t[2] if len(t)>2 and t[2] is not None else t[1] if len(t)>1 and t[1] is not None else np.inf for t in tuples], dtype=float)
        else:
            self.chroms = np.array(chroms, dtype=int) if chroms is not None else np.repeat("", len(starts))
            self.starts = np.array(starts, dtype=int) if starts is not None else np.repeat(0, len(chroms))
            self.ends = np.array(ends) if ends is not None else np.array(starts) if starts is not None else np.repeat(np.inf, len(chroms))
            assert len(self.starts)==len(self.ends)==len(self.chroms)
    
    def containsPoint(self, pos, chrom=""):
        return (self.chroms == chrom) & (self.starts <= pos) & (pos <= self.ends)
    
    def asRegionText(self):
        return ["{}{}{}{}{}".format(self.chroms[i],
                                    ":" if self.chroms[i] != "" and self.starts[i] > 0 else "",
                                    int(self.starts[i]) if self.starts[i] > 0 else "",
                                    "-" if self.starts[i] > 0 and self.ends[i] < np.inf else "",
                                    int(self.ends[i]) if self.starts[i] > 0 and self.ends[i] < np.inf else "") for i in range(len(self.starts))]

