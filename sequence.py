#!/usr/bin/env python

#script that can read and write sequence files in fasta and phylip format
# reads from stdin and writes to stdout

import sys
import argparse
import genomics

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--phylipIn", help="Input is phylip format", action = "store_true")
parser.add_argument("-P", "--phylipOut", help="Output is phylip format", action = "store_true")

parser.add_argument("-r", "--regions",nargs="+", action = "store", metavar="region",
                    help="Output regions and orientation e.g. 'SEQX:1001-1500:+' (FROM, TO and orientation are optional)")

parser.add_argument("-f", "--regionsFile", action = "store",
                    help="File of regions to output. Seq Name and optionally from and to and orintation.Tab separated")

parser.add_argument("-l", "--lineLen", help="Output line length", type = int, action = "store", metavar = "integer", default = 100)

parser.add_argument("--extendLeft", help="Extend left by", type = int, action = "store", metavar = "integer", default=0)
parser.add_argument("--extendRight", help="Extend right by", type = int, action = "store", metavar = "integer", default=0)

parser.add_argument("--truncateNames", help="Truncate names at first whitespace", action = "store_true")

parser.add_argument("--preserveNames", action = "store_true",
                    help = "Do not add start and end position to names of chopped sequences")


args = parser.parse_args()
#args = parser.parse_args("-n 5 -t test.trees -o test.topos.txt -w test.weights.B.csv -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

if args.phylipIn: inFormat = "phylip"
else: inFormat = "fasta"

if args.phylipOut: outFormat = "phylip"
else: outFormat = "fasta"

l = args.lineLen


########################################################################

allText = sys.stdin.read()

if inFormat == "fasta": names, seqs = genomics.parseFasta(allText)
else: names, seqs = genomics.parsePhylip(allText)

if args.truncateNames: names = [name.split()[0] for name in names]

regions = [genomics.parseRegionText(r) for r in args.regions] if args.regions else []

if args.regionsFile:
        with open(args.regionsFile, "r") as rf:
            for line in rf: regions.append(genomics.parseRegionList(line.split()))

#only filter and chop sequences if necessary
if len(regions) >= 1:
    outNames = []
    outSeqs = []
    for seqName,start,end,ori in regions:
        i = names.index(seqName)
        outNames.append(seqName)
        if start != None or end != None or ori == "-":
            seqLen = len(seqs[i])
            if start == None: start = 1
            if end == None: end = seqLen
            start = max(1,start-args.extendLeft)
            end = min(seqLen,end+args.extendRight)
            if ori == "-": outSeqs.append(genomics.revTrans(seqs[i][start-1:end]))
            else: outSeqs.append(seqs[i][start-1:end])
            if not args.preserveNames: outNames[-1] = outNames[-1] + ":" + str(start) + "-" + str(end) + ":" + ori
        else: outSeqs.append(seqs[i])
else:
    outNames = names
    outSeqs = seqs

sys.stderr.write("\nWriting %i sequences.\n" %len(outNames))

sys.stdout.write(genomics.makeAlnString(names=outNames,seqs=outSeqs,outFormat=outFormat,lineLen=l))


