#!/usr/bin/env python

import gzip, argparse, sys

def mafBlockReader(mafFile):
    line=mafFile.readline()
    while line[0] != "a": line=mafFile.readline()
    block = []
    while line != "":
        line=mafFile.readline()
        if line == "" or line[0] == "a":
            yield block
            block = []
        elif line[0] == "s": block.append(line)


def parseMafBlock(mafBlock):
    output = {}
    for line in mafBlock:
        source,start,size,strand,srcSize,seq = line.split()[1:]
        output[source] = {"start":int(start), "size":int(size), "strand":strand, "srcSize":int(srcSize), "seq":seq}
    return output


parser = argparse.ArgumentParser()

parser.add_argument("-m", "--mafFile", help="Input MAK file", action = "store", required = False)
parser.add_argument("-g", "--genoFile", help="Output geno file", action = "store", required = False)
parser.add_argument("--ref", help="Sequences to use as reference", action = "store", type=str, required=True)
parser.add_argument("--renameChromAs", help="Value for CHROM column if not ref", action = "store", type=str, required=False)
parser.add_argument("--seqNames", help="Sequences to output", action = "store", nargs="+", type=str, required=True)
parser.add_argument("--renameSeqsAs", help="If you want to rename the sequence columns", nargs="+", type=str, required=False)
parser.add_argument("--minSeqsRequired", help="Minimum sequences required to output block", action = "store", type=int, default=1)
parser.add_argument("--minSize", help="Minimum block size", action = "store", type=int, default=1)
parser.add_argument("--keepLowercase", help="Do not convert lowercase to uppercase", action = "store_true")
parser.add_argument("--lowercaseToN", help="Convert lowercase to N", action = "store_true")

args = parser.parse_args()

############################# read files

if args.mafFile:
    if args.mafFile[-3:] == ".gz": mafFile = gzip.open(args.mafFile, "rt")
    else: mafFile = open(args.mafFile, "rt")
else: mafFile = sys.stdin

if args.genoFile:
    if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "w")
    else: genoFile = open(args.genoFile, "w")
else: genoFile = sys.stdout

if args.renameSeqsAs:
    assert len(args.renameSeqsAs) == len(args.seqNames), "Incorrect number of new sequence names."
    outputSeqNames = args.renameSeqsAs
else: outputSeqNames = args.seqNames

genoFile.write("#CHROM\tPOS\t" + "\t".join(outputSeqNames) + "\n")

chrom = args.renameChromAs if args.renameChromAs else args.ref

############################# make translations

#translation for conversion of gap to N and cleanup lowercase
if args.keepLowercase: cleanupTrans = str.maketrans("-", "N")
elif args.lowercaseToN: cleanupTrans = str.maketrans("-acgtkmryvhbdn", "NNNNNNNNNNNNNN")
else: cleanupTrans = str.maketrans("-acgtkmryvhbdn", "NACGTKMRYVHBDN")

#translation table for complementation
complementTrans = str.maketrans("AaCcGgTtKkMmRrYyVvHhBbDdNn", "TtGgCcAaMmKkYyRrBbDdVvHhNn")

#############################


blockGen = mafBlockReader(mafFile)

for block in blockGen:
    blockData = parseMafBlock(block)
    
    seqNamesPresent = blockData.keys()
        
    sys.stderr.write("\nProcessing block with {} sequences:\n".format(len(seqNamesPresent)))
    for seqName in seqNamesPresent:
        sys.stderr.write("source={}, start={}, size={}, strand={}\n".format(seqName,
                                                                              blockData[seqName]["start"],
                                                                              blockData[seqName]["size"],
                                                                              blockData[seqName]["strand"]))
    
    if args.ref not in seqNamesPresent:
        sys.stderr.write("Reference absent - skipping block.\n")
        continue
    
    if blockData[args.ref]["size"] < args.minSize:
        sys.stderr.write("Block too short - skipping block.\n")
        continue
    
    #if any sequences are missing, record those as empty
    desiredSeqNamesPresent = [seqName for seqName in seqNamesPresent if seqName in args.seqNames]
    sys.stderr.write("{} of {} desired sequences are present\n".format(len(desiredSeqNamesPresent), len(args.seqNames)))
    if len(desiredSeqNamesPresent) < args.minSeqsRequired:
        sys.stderr.write("Too few sequences - skipping block.\n")
        continue
    
    if len(desiredSeqNamesPresent) < len(args.seqNames):
        empty = ["N"]*blockData[args.ref]["size"]
        for seqName in args.seqNames:
            if seqName not in desiredSeqNamesPresent: sequences[seqName] = empty
    
    #record indices that are not gaps in the reference sequence
    
    refTrueLen = blockData[args.ref]["size"]
    
    #alignment length can be longer if there are gaps
    refAlnLen = len(blockData[args.ref]["seq"])
    
    #get corresponding indices in the complete reference alignment (ie all that are not gaps)
    refIndicesInAln = [i for i in range(refAlnLen) if blockData[args.ref]["seq"][i] != "-"]
    
    sequences = {}
    
    #use the reference sequence to determine the positions
    if blockData[args.ref]["strand"] == "-":
        #if reference is reversed, reverse positions
        positions = range(blockData[args.ref]["start"]+1, blockData[args.ref]["start"]+1 - refTrueLen, -1)[::-1]
        #reverse complement all sequences
        for seqName in desiredSeqNamesPresent:
            sequences[seqName] = blockData[seqName]["seq"].translate(cleanupTrans).translate(complementTrans)[::-1]
    else:
        positions = range(blockData[args.ref]["start"]+1, blockData[args.ref]["start"]+1 + refTrueLen)
        for seqName in desiredSeqNamesPresent:
            sequences[seqName] = blockData[seqName]["seq"].translate(cleanupTrans)
    
    #write sequences, but only for non gaps in reference
    for i in range(refTrueLen):
        genoFile.write("\t".join([chrom, str(positions[i]), "\t".join([sequences[seqName][refIndicesInAln[i]] for seqName in args.seqNames])]) + "\n")

#############################

mafFile.close()
genoFile.close()




