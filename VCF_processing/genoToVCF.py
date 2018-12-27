import sys, gzip, argparse
import genomics
from time import sleep

def makeVCFline(scaffold, position, GTdict, names, refDict=None, genoFormat=None):
    genomeSite = genomics.GenomeSite(genoDict = GTdict, sampleNames=names, genoFormat=genoFormat)
    alleles = genomeSite.alleles(byFreq = True)
    if alleles == []: alleles = ["N"]
    
    if refDict:
        refBase = refDict[scaffold][int(position) - 1]
        if refBase in alleles: alleles.pop(alleles.index(refBase))
        alleles = [refBase] + alleles
    else: refBase = alleles[0]
    
    alt = alleles[1:]
    if alt == []: alt = ["."]
    
    codedGenos = genomeSite.asList(mode="coded", alleles=alleles)
    output = [scaffold,str(position),".",refBase,",".join(alt),".",".",".","GT"] + codedGenos
    return "\t".join(output)


#########################################################################################################################

### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","paired"])
parser.add_argument("-o", "--outFile", help="Output vcf file", action = "store")
parser.add_argument("-r","--reference", help="Reference fasta", action='store')
parser.add_argument("-s", "--samples", help="Analysis threads", action = "store")

args = parser.parse_args()

if args.genoFile:
    genoFile = gzip.open(args.genoFile,"r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin

if args.outFile:
    outFile = gzip.open(args.outFile,"w") if args.outFile.endswith(".gz") else open(args.outFile, "w")
else: outFile = sys.stdout


if args.reference:
    sys.stderr.write("Parsing reference. This could take a while...\n")
    try:
        with open(args.reference + ".fai","r") as fai:
            scafs_lengths = [line.split()[:2] for line in fai]
    except:
        sys.stderr.write("WARNING: Could not parse fai file, vcf header will not contain contig entries...\n")
        scafs_lengths = None
    
    with gzip.open(args.reference,"r") if args.reference.endswith(".gz") else open(args.reference, "r") as ref:
            refDict = dict(zip(*genomics.parseFasta(ref.read())))

else: refDict = None
#########################################################################################

genoFileReader = genomics.GenoFileReader(genoFile)

allNames = genoFileReader.names

if not args.samples: namesToUse = allNames
else: namesToUse = args.samples.split(",")

outFile.write("##fileformat=VCFv4.2\n")

if refDict:
    outFile.write("##reference=file:{}\n".format(args.reference.split("/")[-1]))
    
    if scafs_lengths:
        for x in range(len(scafs_lengths)):
            outFile.write("##contig=<ID={},length={}>\n".format(scafs_lengths[x][0],scafs_lengths[x][1]))

outFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n")

outFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(name for name in namesToUse) + "\n")


linesDone = 0

print >> sys.stderr, "Converting...\n"

for siteData in genoFileReader.siteBySite():
    outFile.write(makeVCFline(siteData["scaffold"], siteData["position"],
                              siteData["GTs"], namesToUse, refDict, args.genoFormat) + "\n")
    linesDone += 1
    if linesDone % 100000 == 0: print >> sys.stderr, linesDone, "lines converted...\n"

genoFile.close()
outFile.close()
