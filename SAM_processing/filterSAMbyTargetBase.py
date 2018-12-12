#!/usr/bin/env python

#script to extract reads from a bam file that have A CERTAIN BASE at A CERTAIN POSITION

import pysam, argparse, sys, gzip, os

parser=argparse.ArgumentParser()

parser.add_argument("-i", "--inBam", help="Input bam file", action = "store", required = True)
parser.add_argument("-o", "--outBam", help="Output bam file", action = "store", required = True)
parser.add_argument("-t", "--targetsFile", help="Targets file (contig, position and base)", action='store', required = True)

#args = parser.parse_args("-i cry.SM16.K570.DSW31521_HF5KHALXX_L8_HF5YNALXX_L2.170110.KIT1701.stampy.s02.rmdup.bam -o cry.SM16.K570.DSW31521_HF5KHALXX_L8_HF5YNALXX_L2.170110.KIT1701.stampy.s02.rmdup.neoW_targets.bam -t ../vcf/queens45.KIT1701.stampy.rmdup.HC.DP8GQ20.fusedPrivateFixed.allele.geno.gz".split())
args = parser.parse_args()


###

#input bam reader
inBam = pysam.AlignmentFile(args.inBam, "rb")

#output bam writer
outBam = pysam.AlignmentFile(args.outBam + "_unsorted", "wb", template = inBam)


#object to store selected read names
selectedReadNames=set()

#read targets file
targetsFile = gzip.open(args.targetsFile, "r") if args.targetsFile.endswith("gz") else open(args.targetsFile, "r")

sys.stderr.write("\nFinding entries containing the target base...\n")
for targetsLine in targetsFile:
    if targetsLine.startswith("#"): continue
    sys.stderr.write(".")
    targetContig,targetPos,targetBase = targetsLine.split()
    targetPos = int(targetPos) - 1
    #fetch reads from this location
    entries = inBam.fetch(contig=targetContig, start=targetPos, stop = targetPos+1)
    for entry in entries:
        try: queryPairedPositions,refPairedPositions = zip(*entry.get_aligned_pairs())
        except: continue
        if targetPos in refPairedPositions:
            readTargetPos = queryPairedPositions[refPairedPositions.index(targetPos)]
            if readTargetPos != None:
                readTargetBase = entry.query_sequence[readTargetPos]
                if readTargetBase.upper() == targetBase:
                    selectedReadNames.add(entry.query_name)


sys.stderr.write("\nFound {} entries carrying a target base\n".format(len(selectedReadNames)))


#index file for pulling out reads
sys.stderr.write("\nBuilding read name index for extracting selected pairs...\n")
readIndex = pysam.IndexedReads(inBam)
readIndex.build()

sys.stderr.write("\nWriting {} selected entries...\n".format(len(selectedReadNames)))
for readName in selectedReadNames:
    entries=readIndex.find(readName)
    for entry in entries: outBam.write(entry)

inBam.close()
outBam.close()
targetsFile.close()

sys.stderr.write("\nSorting output bam file...\n")

pysam.sort("-o", args.outBam, args.outBam + "_unsorted")

os.remove(args.outBam + "_unsorted")

sys.stderr.write("\nIndexing sorted bam file...\n")

pysam.index(args.outBam)

sys.stderr.write("\nDone.\n")

