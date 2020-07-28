#!/usr/bin/env python

#20 July 2017
#Originaly written to oconvert gff positions,
#now it works on any generic genomic file with scaffold and position (and optionally end position)
#This goes hand-in-hand with transferFasta.py, which will rebuild a fasta references from an old one
#agp option is untested, but should work the same as providing a "transfers" file

import sys, argparse, gzip, time
import numpy as np

def contains(array, x):
    return (x >= array[:,0]) & (x <= array[:,1])

def newPos(pos, start=1, newStart=None, newEnd=None, reverse=False,):
    pos = pos-start+1
    if not reverse:
        assert newStart is not None, "newStart must be provided."
        return newStart+pos-1
    else:
        assert newEnd is not None, "newEnd must be provided."
        return  newEnd-pos+1

##################################################

parser=argparse.ArgumentParser()

parser.add_argument("-i", "--inFile", help="Input file", action = "store")
parser.add_argument("-o", "--outFile", help="Output gff file", action = "store")
parser.add_argument("-p", "--preset", help="Preset file type", action = "store", choices = ("vcf","gff",))
parser.add_argument("--scafCol", help="Column of scaffold name", action = "store", type=int, default=1)
parser.add_argument("--startCol", help="Column of start position", action = "store", type=int, default=2)
parser.add_argument("--endCol", help="Column of end position", action = "store", type=int, default=2)
parser.add_argument("--strandCol", help="Column of strand (orientation)", action = "store", type=int)
parser.add_argument("--sep", help="Input file separator", action = "store", default=None)
parser.add_argument("-f", "--failFile", help="Failed lines file", action = "store")
parser.add_argument("-a", "--agpFile", help="AGP file for position conversion", action='store')
parser.add_argument("-t", "--transfersFile", help="Chrom and position transfer table (alternative to AGP)", action='store')
parser.add_argument("--header", help="File has a header line", action = "store_true")
parser.add_argument("--keepFails", help="Retain failed transfers in output file", action = "store_true")
parser.add_argument("--allowAGPfails", help="Skip malformed lines in the agp without quitting", action = "store_true")

args = parser.parse_args()

if not args.inFile: inFile = sys.stdin
else: inFile = gzip.open(args.inFile,"rt") if args.inFile.endswith(".gz") else open(args.inFile,"rt")

if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"wt") if args.outFile.endswith(".gz") else open(args.outFile,"wt")

if args.failFile:
    failFile = gzip.open(args.failFile,"wt") if args.failFile.endswith(".gz") else open(args.failFile,"wt")
else:
    failFile = open('/dev/null', 'wt')
    if not args.keepFails:
        sys.stderr.write("\nWARNING: Failed transfers will not be shown. To catch them, specify a --failFile\n\n") 


if args.preset == "vcf": scafCol, startCol, endCol, strandCol = (1,2,2,None)
elif args.preset == "gff":scafCol, startCol, endCol, strandCol = (1,4,5,7)
else:scafCol, startCol, endCol, strandCol = (args.scafCol, args.startCol, args.endCol, args.strandCol)


sep=args.sep
outsep = sep if sep != None else "\t"

if not args.transfersFile and not args.agpFile:
    raise ValueError("Please provide an AGP file (or a 'transfers' file)")

###################################################
#get transfer info

transfers = {}

if args.agpFile:
    with open(args.agpFile, "rt") as agpFile:
        for line in agpFile:
            if not line.startswith("#"):
                try:newScaf,newStart,newEnd,part,component,scaf,start,end,strand = line.split()
                except:
                    if args.allowAGPfails:
                        sys.stderr.write("WARNING: skipping malformed agp line:\n{}".format(line))
                        continue
                    else:
                        raise ValueError("agp file should have nine fields.")
                if component == "N" or component == "U": continue
                if scaf not in transfers: transfers[scaf] = np.empty(shape = [0,3])
                transfers[scaf] = np.vstack([transfers[scaf], np.array([(int(start), int(end),
                                                                        {"scaf":scaf, "start":int(start), "end":int(end),"strand":strand,
                                                                        "newScaf":newScaf,"newStart":int(newStart),"newEnd":int(newEnd)})])])

else:
    with open(args.transfersFile, "rt") as transfersFile:
        for line in transfersFile:
            if not line.startswith("#"):
                try:newScaf,newStart,newEnd,scaf,start,end,strand = line.split()
                except:
                    raise ValueError("Transfers file should have seven fields for newChrom, newStart, newEnd, chrom, start, end and strand.")
                if scaf not in transfers: transfers[scaf] = np.empty(shape = [0,3])
                transfers[scaf] = np.vstack([transfers[scaf], np.array([(int(start), int(end),
                                                                        {"scaf":scaf, "start":int(start), "end":int(end),"strand":strand,
                                                                        "newScaf":newScaf,"newStart":int(newStart),"newEnd":int(newEnd)})])])

###################################################

#transfer chrom and positions

if args.header:
    headLine = inFile.readline()
    outFile.write(headLine)
    failFile.write(headLine)

for line in inFile:
    if not line.startswith("#"):
        elements = line.strip().split(sep)
        scaf = elements[scafCol-1]
        start = int(elements[startCol-1])
        end = int(elements[endCol-1])
        strand = elements[strandCol-1] if strandCol else "+"
        assert strand is "+" or strand is "-"
        if scaf in transfers:
            startInterval = np.where(contains(transfers[scaf], start))[0]
            endInterval = np.where(contains(transfers[scaf], end))[0]
            if len(startInterval) == len(endInterval) == 1 and startInterval[0] == endInterval[0]:
                intervalData = transfers[scaf][startInterval[0],:][2]
                newScaf = intervalData["newScaf"]
                if intervalData["strand"] == "+":
                    newStart = str(newPos(pos=start, start=intervalData["start"],
                                              newStart=intervalData["newStart"], newEnd=intervalData["newEnd"], reverse=False))
                    newEnd = str(newPos(pos=end, start=intervalData["start"],
                                            newStart=intervalData["newStart"], newEnd=intervalData["newEnd"], reverse=False))
                    newStrand = strand
                else:
                    newStart = str(newPos(pos=end, start=intervalData["start"],
                                              newStart=intervalData["newStart"],newEnd=intervalData["newEnd"], reverse=True))
                    newEnd = str(newPos(pos=start, start=intervalData["start"],
                                            newStart=intervalData["newStart"], newEnd=intervalData["newEnd"], reverse=True))
                    newStrand = "-" if strand == "+" else "+"
                
                elements[scafCol-1] = newScaf
                elements[startCol-1] = newStart
                elements[endCol-1] = newEnd
                if strandCol: elements[strandCol-1] = newStrand
                outFile.write(outsep.join(elements) + "\n")
            
            else:
                failFile.write("#BROKEN\n")
                failFile.write(outsep.join(elements) + "\n")
                
                if args.keepFails:
                    elements[scafCol-1] = "NA"
                    elements[startCol-1] = "NA"
                    elements[endCol-1] = "NA"
                    if strandCol: elements[strandCol-1] = "NA"
                    outFile.write(outsep.join(elements) + "\n")
                
        else:
            failFile.write("#MISSING\n")
            failFile.write(outsep.join(elements) + "\n")
            
            if args.keepFails:
                elements[scafCol-1] = "NA"
                elements[startCol-1] = "NA"
                elements[endCol-1] = "NA"
                if strandCol: elements[strandCol-1] = "NA"
                outFile.write(outsep.join(elements) + "\n")

outFile.close()
failFile.close()

