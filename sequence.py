#!/usr/bin/python

#script that can read and write sequence files in fasta and phylip format
# reads from stdin and writes to stdout

import sys
import argparse
import DNA

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--phylipIn", help="Input is phylip format", action = "store_true")
parser.add_argument("-P", "--phylipOut", help="Output is phylip format", action = "store_true")

parser.add_argument("-n", "--seqNames", help="Output only these sequences (comma separated)", action = "store", metavar = "names")

parser.add_argument("-l", "--lineLen", help="Output line length", type = int, action = "store", metavar = "integer", default = 100)

parser.add_argument("--truncateNames", help="Truncate names at first whitespace", action = "store_true")


args = parser.parse_args()
#args = parser.parse_args("-n 5 -t test.trees -o test.topos.txt -w test.weights.B.csv -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

if args.phylipIn: inFormat = "phylip"
else: inFormat = "fasta"

if args.phylipOut: outFormat = "phylip"
else: outFormat = "fasta"

l = args.lineLen

if args.seqNames: outNames = args.seqNames.split(",")

########################################################################

allText = sys.stdin.read()

if inFormat == "fasta": names, seqs = DNA.parseFasta(allText)
else: names, seqs = DNA.parsePhylip(allText)

if args.truncateNames: names = [name.split()[0] for name in names]

if outNames:
    keep = [names.index(name) for name in outNames]
    names = [names[k] for k in keep]
    seqs = [seqs[k] for k in keep]
    
sys.stdout.write(DNA.makeAlnString(names,seqs,outFormat,l))


