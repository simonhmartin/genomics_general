#!/usr/bin/env python

import sys,argparse,gzip
import numpy as np
import genomics
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", action = "store",
                    choices = ["haplo", "diplo", "pairs", "alleles", "phased"], default="phased",
                    help="Genotype format [otherwise will be inferred (slower)]")
parser.add_argument("--prefix", help="Output file prefix for PED and MAP files", action = "store")
parser.add_argument("--makeFAM", help="Make FAM file (assumes all unrelated)", action = "store_true")
parser.add_argument("--FAMprefix", help="Output file prefix for FAM file", action = "store")
parser.add_argument("-s", "--samples", help="sample names", nargs="+", action='store')


args = parser.parse_args()

#args = parser.parse_args(["-g", "test.geno", "-f", "phased"])


if args.genoFile:
    genoFile = gzip.open(args.genoFile,"r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else:
    assert args.prefix != None, "Please povide a prefix for the ouput files"
    genoFile = sys.stdin


#########################################################################################
###################### read data ####################################################


#we will make a list of geno windows for each scaffold
scafWindows = []
for i, window in enumerate(genomics.nonOverlappingSitesWindows(genoFile, windSites=np.inf,names=args.samples)):
    scafWindows.append(window)
    sys.stderr.write("{} scaffolds read into memory\n".format(i+1))

genoFile.close()


#make concatenated sequences for each sample
names = scafWindows[0].names
sampleSeqs = defaultdict(list)
for scafWindow in scafWindows:
    windowSeqs = scafWindow.seqDict()
    for name in names:
        sampleSeqs[name] += [a for pair in zip(*genomics.splitSeq(windowSeqs[name],genoFormat=args.genoFormat)) for a in pair]


#########################################################################################
###################### write outputs ####################################################

prefix = args.prefix if args.prefix else args.genoFile.rsplit(".", 1)[0]

sys.stderr.write("Writing PED file...\n")

#write sequences for each sample
with open(prefix + ".ped" , "w") as outPed:
    for name in scafWindows[0].names:
        outPed.write(" ".join(["0", name, "0 0 0 0 "]))
        outPed.write(" ".join(sampleSeqs[name]) + "\n")

sys.stderr.write("Writing MAP file...\n")

#write map file
with open(prefix + ".map" , "w") as outMap:
    for scafWindow in scafWindows:
        for pos in scafWindow.positions:
            outMap.write("{} {} 0 {}\n".format(scafWindow.scaffold, pos, pos))

if args.makeFAM:
    sys.stderr.write("Writing FAM file...\n")
    #write fam file
    with open(args.FAMprefix if args.FAMprefix else prefix + ".fam" , "w") as outFam:
        for name in names:
            outFam.write("0 {} 0 0 0 0\n".format(name))
