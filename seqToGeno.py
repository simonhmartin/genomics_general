#!/usr/bin/env python

import gzip, argparse, sys

import genomics

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--seqFile", help="Input sequence file", action = "store", required = False)
parser.add_argument("-g", "--genoFile", help="Output geno file", action = "store", required = False)
parser.add_argument("-f", "--format", help="Sequence file format", action = "store", required = False, choices = ("phylip","fasta",), default = "fasta")
parser.add_argument("-M", "--mode", help="Output mode for separate sequences", action = "store", choices = ("samples","contigs"), required = False, default = "samples")
parser.add_argument("-C", "--chrom", help="Name for chromosome or contig if sequences are samples", action = "store", required = False, default = "contig0")
parser.add_argument("-N", "--name", help="Name for sample if sequences are contigs", action = "store", required = False, default = "sample0")
parser.add_argument("-S", "--sequences", help="Sequences to output", action = "store", nargs="+", type=str)
parser.add_argument("--merge", help="For multi-phylip input, do not increment scaffold numbers", action = "store_true")
parser.add_argument("-P", "--ploidy", help="Ploidy for joining sequences", action = "store", nargs="+", type=int, default = [1])
parser.add_argument("--randomPhase", help="Randomize phase for fused sequences", action = "store_true")

args = parser.parse_args()

if args.seqFile:
    if args.seqFile[-3:] == ".gz": seqFile = gzip.open(args.seqFile, "rt")
    else: seqFile = open(args.seqFile, "r")
else: seqFile = sys.stdin

if args.genoFile:
    if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "wt")
    else: genoFile = open(args.genoFile, "w")
else: genoFile = sys.stdout

if args.randomPhase:
    import random

#############################

#read sequence file
seqString = seqFile.read()

#parse
if args.format == "fasta":
    seqNames,seqs = genomics.parseFasta(seqString)
    multi = False

elif args.format == "phylip":
    #with phylip its possible to have multiple alignments, so we need to check if thats the case
    pieces = genomics.parsePhylip(seqString)
    if type(pieces) == tuple:
        seqNames,seqs = pieces
        multi = False
    else:
        _seqNames_,_seqs_ = zip(*pieces)
        multi = True


if not multi:
    #if there is a single set of sequences we parse it and output either as contigs or individuals
    #sequences to keep
    if args.sequences is not None:
        seqs = [seqs[seqNames.index(x)] for x in args.sequences]
        seqNames = args.sequences

    #phase if ploidy is not 1
    if max(args.ploidy) > 1:
        seqs, seqNames = genomics.haploToPhased(seqs, seqNames=seqNames, ploidy = args.ploidy, randomPhase=args.randomPhase)

    ###output
    
    if args.mode == "samples":
        genoFile.write("#CHROM\tPOS\t" + "\t".join(seqNames) + "\n")
        for x in range(len(seqs[0])): genoFile.write(args.chrom + "\t" + str(x+1) + "\t" + "\t".join([s[x] for s in seqs]) + "\n")
    
    elif args.mode == "contigs":
        genoFile.write("#CHROM\tPOS\t" + args.name + "\n")
        for y in range(len(seqNames)):
            for x in range(len(seqs[y])): genoFile.write(seqNames[y] + "\t" + str(x+1) + "\t" + seqs[y][x] + "\n")

else:
    #otherwise we output each set as a contig, with each sequence as a sample
    #first confirm all same length
    assert len(set(map(len,_seqNames_))) == 1, "For multi phylip, all alignments must have same number of sequences"
    #sequences to keep
    seqNames = args.sequences if args.sequences else _seqNames_[0]
    #now we need to stich together all sequences, but they are possibly different orders, so index each one each time
    _indices_ = [[_names_.index(name) for name in seqNames] for _names_ in _seqNames_]
    _seqs_ = [[_seqs_[i][j] for j in _indices_[i]] for i in range(len(_seqs_))]
    
    #phase if ploidy is not 1
    if max(args.ploidy) > 1:
        _seqs_ = [genomics.haploToPhased(seqs, ploidy = args.ploidy) for seqs in _seqs_]
        seqNames = genomics.makePhasedNames(seqNames, ploidy = args.ploidy, randomPhase=args.randomPhase)
    
    ###output
    #with multi data there is only one output mode, so we ignore the mode specified.
    genoFile.write("#CHROM\tPOS\t" + "\t".join(seqNames) + "\n")
    for i,seqs in enumerate(_seqs_):
        contigName = args.chrom if args.merge else args.chrom + str(i)
        for x in range(len(seqs[0])): genoFile.write(contigName + "\t" + str(x+1) + "\t" + "\t".join([s[x] for s in seqs]) + "\n")

#############################


seqFile.close()
genoFile.close()
