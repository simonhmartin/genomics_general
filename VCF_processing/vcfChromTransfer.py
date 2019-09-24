#!/usr/bin/env python

import sys, argparse, gzip, subprocess, string
from collections import OrderedDict

def newPos(pos, start=1, newStart=None, newEnd=None, reverse=False,):
    pos = pos-start+1
    if not reverse:
        assert newStart is not None, "newStart must be provided."
        return newStart+pos-1
    else:
        assert newEnd is not None, "newEnd must be provided."
        return  newEnd-pos+1


#improved tabix stream that fixes the bug where readining is terminated prematurely
def tabixStream(fileName, region = None, chrom = None, start=None, end=None, header=False, iterable=True):
    if not region: region = chrom+":"+str(start)+"-"+str(end)  
    
    if header:
        p = subprocess.Popen(['tabix', '-h', fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    else:
        p = subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    
    if iterable: return iter(p.communicate()[0].strip().split("\n"))
    else: return p.communicate()[0].strip().split("\n")


#translation table for bases
if sys.version_info>=(3,0):
    complementTrans = str.maketrans("ACGT", "TGCA")
else:
    complementTrans = string.maketrans("ACGT", "TGCA")

##################################################

parser=argparse.ArgumentParser()

parser.add_argument("-v", "--vcfFile", help="Input vcf file (tabix indexed)", action = "store", required = True)
parser.add_argument("-o", "--outFile", help="Output vcf file", action = "store")

parser.add_argument("-a", "--agpFile", help="AGP file for position conversion", action='store')
parser.add_argument("-t", "--transfersFile", help="Chromosome and position transfer table", action='store')

parser.add_argument("--chroms", help="If only specific output chromosomes are desired", nargs="+", action='store')

args = parser.parse_args()

if not args.transfersFile and not args.agpFile:
    raise ValueError("Please provide an AGP file (or a 'transfers' file)")

if not args.outFile: outStream = sys.stdout
else: outStream = gzip.open(args.outFile,"r") if args.outFile.endswith(".gz") else open(args.outFile,"r")

###################################################
#get transfer info

transfers = []

if args.agpFile:
    with open(args.agpFile, "r") as agpFile:
        for line in agpFile:
            if not line.startswith("#"):
                try:newChrom,newStart,newEnd,part,component,chrom,start,end,strand = line.split()
                except:
                    sys.stderr.write("\nWARNING: failed to extract 9 fields from agp line\n{}Line will be ignored\n.".format(line))
                    continue
                if component == "N" or component == "U": continue
                if not args.chroms or newChrom in args.chroms:
                    transfers.append([newChrom,newStart,newEnd,chrom,start,end,strand])

else:
    with open(args.transfersFile, "r") as transfersFile:
        for line in transfersFile:
            if not line.startswith("#"):
                try:newChrom,newStart,newEnd,chrom,start,end,strand = line.split()
                except:
                    sys.stderr.write("\nWARNING: failed to extract 7 fields from transfers line\n{}Line will be ignored\n.".format(line))
                    continue
                if not args.chroms or newChrom in args.chroms:
                    transfers.append([newChrom,newStart,newEnd,chrom,start,end,strand])

#the end of the last transfer interval for each chromosome is taken as the chrom length
newChromLengths = OrderedDict()
for t in transfers:
    end = int(t[2])
    if t[0] in newChromLengths:
        assert end > newChromLengths[t[0]], "Transfers for chrom {} not in correct order. {} is not more than than {}\n".format(t[0], end, newChromLengths[t[0]])
    newChromLengths[t[0]] = end

newChroms = newChromLengths.keys()

###################################################
#read and write vcf header, replacing the chrom info
contigs_output = False

with gzip.open(args.vcfFile, "rt") as vcf:
    for line in vcf:
        if line.startswith('##contig'):
            if not contigs_output:
                for c in newChroms: outStream.write('##contig=<ID={},length={}>\n'.format(c, str(newChromLengths[c])))
                contigs_output = True
            continue
        elif line.startswith('#'):
            outStream.write(line)
            continue
        else: break

###################################################

#extract each interval from the vcf and transfer chrom and positions

for t in transfers:
    newChrom,newStart,newEnd,chrom,start,end,strand = t
    region = chrom+":"+start+"-"+end
    sys.stderr.write("\nGetting region " + region + " from vcf...\n") 
    vcfLines = tabixStream(args.vcfFile, region, iterable=False) #this might take a while if it's a big chunk
    
    if vcfLines == ['']:
        sys.stderr.write("WARNING: Region empty. If this is unexpected, ensure input vcf is tabix indexed.\n") 
        continue
    
    sys.stderr.write("Region extracted. {} lines.\n".format(len(vcfLines))) 
    
    if strand == "-":
        reverse = True
        sys.stderr.write("Orientation is reverse.\n")
        sys.stderr.write("reversing...\n")
        vcfLines.reverse()
    else:
        reverse = False
        sys.stderr.write("Orientation is forward.\n")
    
    sys.stderr.write("Writing new region " + newChrom+":"+newStart+"-"+newEnd + "...\n") 
    for vcfLine in vcfLines:
        fields = vcfLine.split("\t")
        assert fields[0] == chrom, "Something went wrong: Found chrom {} but expected chrom {}.".format(fields[0],chrom)
        fields[0] = newChrom
        #get new position
        fields[1] = str(newPos(int(fields[1]), start=int(start), newStart=int(newStart), newEnd=int(newEnd), reverse=reverse))
        #complement ref and alternate alleles
        fields[3] = fields[3].translate(complementTrans)
        fields[4] = fields[4].translate(complementTrans)
        outStream.write("\t".join(fields)+"\n")

outStream.close()
