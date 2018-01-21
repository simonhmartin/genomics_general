import gzip, argparse, sys


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

def chunkList(l, nChunks = None, chunkSize = None, return_indices=False):
    N = len(l)
    assert not nChunks is chunkSize is None
    if nChunks is not None:
        assert N % nChunks == 0, "list must be divizable by number of chunks"
        chunkSize = [N/nChunks]*nChunks
    elif isinstance(chunkSize, int):
        assert N % chunkSize == 0, "list must be divizable by chunk size"
        chunkSize = [chunkSize]*(N/chunkSize)
    else: assert N == sum(chunkSize), "Chunk sizes must sum to list length"
        
    indices = []
    r = range(N)
    i = 0
    for c in chunkSize:
        indices.append(range(i,i+c))
        i = i+c
    if return_indices: return ([[l[x] for x in ind] for ind in indices], indices,)
    else: return [[l[x] for x in ind] for ind in indices]



#############################

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--seqFile", help="Input sequence file", action = "store", required = False)
parser.add_argument("-g", "--genoFile", help="Output geno file", action = "store", required = False)
parser.add_argument("-f", "--format", help="Sequence file format", action = "store", required = False, choices = ("phylip","fasta",), default = "fasta")
parser.add_argument("-M", "--mode", help="Output mode for separate sequences", action = "store", choices = ("samples","contigs"), required = False, default = "samples")
parser.add_argument("-C", "--chrom", help="Name for chromosome or contig if sequences are samples", action = "store", required = False, default = "contig0")
parser.add_argument("-N", "--name", help="Name for sample if sequences are contigs", action = "store", required = False, default = "sample0")
parser.add_argument("-S", "--sequences", help="Sequences to output", action = "store", nargs="+", type=str)
parser.add_argument("-P", "--ploidy", help="Ploidy for joining sequences", action = "store", nargs="+", type=int, default = 1)
parser.add_argument("--randomPhase", help="Randomize phase for fused sequences", action = "store_true")

args = parser.parse_args()

if args.seqFile:
    if args.seqFile[-3:] == ".gz": seqFile = gzip.open(args.seqFile, "r")
    else: seqFile = open(args.seqFile, "r")
else: seqFile = sys.stdin

if args.genoFile:
    if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "w")
    else: genoFile = open(args.genoFile, "w")
else: genoFile = sys.stdout

if args.randomPhase:
    import random

#############################

#read sequence file
seqString = seqFile.read()

#parse
if args.format == "fasta": seqNames,seqs = parseFasta(seqString)
elif args.format == "phylip": seqNames,seqs = parsePhylip(seqString)

#sequences to keep

if args.sequences is not None:
    seqs = [seqs[seqNames.index(x)] for x in args.sequences]
    seqNames = args.sequences


#phase if ploidy is not 1
if args.ploidy is not 1:
    ploidy = args.ploidy
    if len(ploidy) == 1:
        assert len(seqNames) % ploidy[0] == 0, "Sequence number must be divizable by ploidy"
        ploidy = ploidy*(len(seqNames)/ploidy[0])
    else: assert len(seqNames) == sum(ploidy), "Ploidys must sum to number of sequences"
    
    indices = chunkList(range(len(seqNames)), chunkSize=ploidy[0], return_indices=True)[1]
    
    zipSeqs = [zip(*[seqs[x] for x in ind]) for ind in indices]
    #randomize phase if necessary
    if args.randomPhase:
        for i in range(len(indices)):
            if ploidy[i] > 1:
                for j in range(len(zipSeqs[i])):
                    zipSeqs[i][j] = random.sample(zipSeqs[i][j], ploidy[i])
    seqs = [["|".join(x) for x in zipSeq] for zipSeq in zipSeqs]
    seqNames = ["_".join([seqNames[x] for x in ind]) for ind in indices]


#############################


if args.mode == "samples":
    genoFile.write("#CHROM\tPOS\t" + "\t".join(seqNames) + "\n")
    for x in xrange(len(seqs[0])): genoFile.write(args.chrom + "\t" + str(x+1) + "\t" + "\t".join([s[x] for s in seqs]) + "\n")

elif args.mode == "contigs":
    genoFile.write("#CHROM\tPOS\t" + args.name + "\n")
    for y in range(len(seqNames)):
        for x in xrange(len(seqs[y])): genoFile.write(seqNames[y] + "\t" + str(x+1) + "\t" + seqs[y][x] + "\n")
    

#############################


seqFile.close()
genoFile.close()

