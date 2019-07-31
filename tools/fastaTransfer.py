#!/usr/bin/env python

#Script to construct new scaffolds from chunks of old ones
#This goes hand-in-hand with transferScafPos.py, which will convert column files with genome coordianates

import sys, argparse, gzip, string
import numpy as np

def parseFasta(string):
    splitString = string.split(">")[1:]
    names = [s.split()[0] for s in splitString]
    seqs = [s[s.index("\n"):].replace("\n","").replace(" ","") for s in splitString]
    return (names,seqs)

seqtrans = string.maketrans("ACGTKMRYNacgtkmryn", "TGCAMKYRNtgcamkyrn")

def revTrans(seq):
    return seq.translate(seqtrans)[::-1]

def subset(things,subLen):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]

def makeFastaString(names=None, seqs=None, seqDict=None, lineLen=60):
    if seqDict: names, seqs = zip(*seqDict.items())
    else: assert len(names) == len(seqs)
    seqs = ["".join(s) for s in seqs]
    output = []
    nSamp = len(names)
    seqLen = max(map(len,seqs))
    if lineLen: seqs = ["\n".join(subset(s,lineLen)) for s in seqs]
    for x in range(nSamp):
        output.append(">" + names[x])
        output.append(seqs[x])
    
    return "\n".join(output) + "\n"

##################################################

parser=argparse.ArgumentParser()

parser.add_argument("-i", "--inFile", help="Input fasta file", action = "store")
parser.add_argument("-o", "--outFile", help="Output fasta file", action = "store")
parser.add_argument("-a", "--agpFile", help="AGP file for position conversion", action='store')
parser.add_argument("-t", "--transfersFile", help="Chromosome and position transfer table", action='store')

args = parser.parse_args()

if not args.inFile: inFile = sys.stdin
else: inFile = gzip.open(args.inFile,"r") if args.inFile.endswith(".gz") else open(args.inFile,"r")

if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"w") if args.outFile.endswith(".gz") else open(args.outFile,"w")

if not args.transfersFile and not args.agpFile:
    raise ValueError("Please provide an AGP file (or a 'transfers' file)")

###################################################

#record all pieces for each new scaffold

pieces = {}

newScafs = []

if args.agpFile:
    with open(args.transfersFile, "r") as transfersFile:
        for line in transfersFile:
            if not line.startswith("#"):
                try:newScaf,newStart,newEnd,part,component,scaf,start,end,strand = line.split()
                except:
                    raise ValueError("AGP file should have 9 columns")
                if component == "N" or component == "U": continue
                if newScaf not in newScafs:
                    newScafs.append(newScaf)
                    pieces[newScaf] = np.empty(shape = [0,3])
                try: pieces[newScaf] = np.vstack([pieces[newScaf], np.array([(int(newStart), int(newEnd),
                                                                            {"scaf":scaf, "start":int(start), "end":int(end),"strand":strand,
                                                                            "newScaf":newScaf,"newStart":int(newStart),"newEnd":int(newEnd)})])])
                except: pass

else:
    with open(args.transfersFile, "r") as transfersFile:
        for line in transfersFile:
            if not line.startswith("#"):
                try:newScaf,newStart,newEnd,scaf,start,end,strand = line.split()
                except:
                    raise ValueError("Transfers file should have seven fields for newChrom, newStart, newEnd, chrom, start, end and strand.")
                if newScaf not in newScafs:
                    newScafs.append(newScaf)
                    pieces[newScaf] = np.empty(shape = [0,3])
                try: pieces[newScaf] = np.vstack([pieces[newScaf], np.array([(int(newStart), int(newEnd),
                                                                            {"scaf":scaf, "start":int(start), "end":int(end),"strand":strand,
                                                                            "newScaf":newScaf,"newStart":int(newStart),"newEnd":int(newEnd)})])])
                except: pass

sys.stderr.write("{} new scaffolds to be made.\n".format(len(newScafs)))


#read fasta and make new seqs
scafs,seqs = parseFasta(inFile.read())
seqDict = dict(zip(scafs,seqs))

newSeqs = []

for newScaf in newScafs:
    sys.stderr.write("Making new sequence: {}, {} pieces, {} bp.\n".format(newScaf,
                                                                         pieces[newScaf].shape[0],
                                                                         np.max(pieces[newScaf][:,1])))
    newSeq = ["N"]*np.max(pieces[newScaf][:,1])
    for piece in pieces[newScaf]:
        pieceData = piece[2]
        pieceSeq = seqDict[pieceData["scaf"]][pieceData["start"]-1:pieceData["end"]]
        if pieceData["strand"] == "-": pieceSeq = revTrans(pieceSeq)
        newSeq[pieceData["newStart"]-1:pieceData["newEnd"]] = pieceSeq
    newSeqs.append(newSeq)


outFile.write(makeFastaString(names=newScafs, seqs = newSeqs))
