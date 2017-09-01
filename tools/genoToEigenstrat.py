import sys,argparse,gzip
import genomics

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","paired"])
parser.add_argument("--genoOutFile", help="Output Eigenstrat geno file", action = "store", required = True)
parser.add_argument("--snpOutFile", help="Output Eigenstrat snp file", action = "store", required = True)
parser.add_argument("-s", "--samples", help="sample names (separated by commas)", action='store')
parser.add_argument("--chromFile", help="Input chromosome number file", action = "store")
parser.add_argument("--cumulativePos", help="For chroms with multiple scafs, add positions to previous scaf", action = "store_true")
parser.add_argument("--nullChrom", help="Chrom_for_unknown_scafs", action = "store", type=int, default=22)


args = parser.parse_args()

if args.genoFile:
    genoFile = gzip.open(args.genoFile,"r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin

genoOut = open(args.genoOutFile, "w")
snpOut = open(args.snpOutFile, "w")

headers = genoFile.readline().split()
allNames = headers[2:]

if args.samples is None: samples = allNames
else: samples = args.samples.split(",")


if args.chromFile:
    chromsProvided = True
    with open(args.chromFile, "r") as chromFile: chromDict = dict([line.split() for line in chromFile]
else: chromDict = {}


linesDone = 0
scaf = None
chrom = None
pos = 0

if args.cumulativePos:
    #dict giving the last known position of each chrom from the previous scaffold
    chromOffset = dict(zip(chromDict.keys(), [0]*len(chromDict)))
    chromOffset[str(args.nullChrom)] = 0


for line in genoFile:
    site = genomics.parseGenoLine(line)
    
    genomeSite = genomics.GenomeSite(genotypes = site.GTs, sampleNames = allNames, genoFormat=args.genoFormat)
    
    if len(genomeSite.alleles()) == 2:
        counts = genomeSite.asList(mode="count", samples=samples, missing = 9)

        genoOut.write("".join([str(c) for c in counts]) + "\n")
        
        if site.scaffold != scaf:
            #different scaffold from the last site
            #if using cumulative positions, change the offset for the last chrom
            if chrom is not None and args.cumulativePos: chromOffset[chrom] = pos
            #now get new scaf, chrom and pos
            scaf = site.scaffold
            try: chrom = chromDict[scaf]
            except: chrom = str(args.nullChrom)
        
        pos = site.position if not args.cumulativePos else site.position + chromOffset[chrom]
        
        snpOut.write(str(linesDone) + "\t" + chrom + "\t" + "0.0" + "\t" + str(pos) + "\n")
    
    linesDone += 1
    if linesDone % 100000 == 0: print linesDone, "lines done..."
    

genoOut.close
snpOut.close()
