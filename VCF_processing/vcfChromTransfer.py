
import sys, argparse, gzip, subprocess
from collections import OrderedDict

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

parser.add_argument("-v", "--vcfFile", help="Input vcf file", action = "store", required = True)
parser.add_argument("-o", "--outFile", help="Output vcf file", action = "store")

parser.add_argument("-t", "--transfersFile", help="Chromosome and position transfer table", action='store', required = True)

args = parser.parse_args()

if not args.outFile: outStream = sys.stdout
else: outStream = gzip.open(args.outFile,"r") if args.outFile.endswith(".gz") else open(args.outFile,"r")

###################################################
#get transfer info

transfersFile = open(args.transfersFile, "r")
assert len(transfersFile.readline().split()) == 7, "Transfers file should have seven fields for newChrom, newStart, newEnd, chrom, start, end and orientation."
transfers = [line.split() for line in transfersFile]
transfersFile.close()

#the end of the last transfer interval is taken as the chrom length
newChromLengths = OrderedDict()
for t in transfers:
    end = int(t[2])
    if t[0] in newChromLengths:
        assert end > newChromLengths[t[0]], "Transfers for each chrom must be ordered such that the last one ends at the chrom end."
    newChromLengths[t[0]] = end

newChroms = newChromLengths.keys()

###################################################
#read and write vcf header, replacing the chrom info
contigs_output = False

vcf = gzip.open(args.vcfFile) if args.vcfFile.endswith(".gz") else open(args.vcfFile)
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

vcf.close()

###################################################

#extract each interval from the vcf and transfer chrom and positions

for t in transfers:
    newChrom,newStart,newEnd,chrom,start,end,ori = t
    region = chrom+":"+start+"-"+end
    sys.stderr.write("Getting region " + region + " from vcf...\n") 
    vcfStream = subprocess.Popen(['tabix',args.vcfFile, region], stdout=subprocess.PIPE)
    vcfLines = vcfStream.stdout.readlines() #this might take a while if it's a big chunk
    
    if ori == "-":
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
        assert fields[0] == chrom, "Something went wrong - tabix pulled out the wrong region."
        fields[0] = newChrom
        fields[1] = str(newPos(int(fields[1]), start=int(start), newStart=int(newStart), newEnd=int(newEnd), reverse=reverse))
        outStream.write("\t".join(fields))

outStream.close()
