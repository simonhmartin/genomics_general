import sys,argparse,gzip,time
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

reader=genomics.GenoFileReader(genoFile)

if args.samples is None: samples = reader.names
else: samples = args.samples.split(",")


if args.chromFile:
    chromsProvided = True
    with open(args.chromFile, "r") as chromFile: chromDict = dict([line.split() for line in chromFile])
else: chromDict = {}


linesDone = 0
scaf = None
chrom = None
pos = 0

if args.cumulativePos:
    #dict giving the last known position of each chrom from the previous scaffold
    chromOffset = dict(zip(chromDict.keys(), [0]*len(chromDict)))
    chromOffset[str(args.nullChrom)] = 0


for siteData in reader.siteBySite():
    
    genomeSite = genomics.GenomeSite(genoDict = siteData["GTs"], genoFormat=args.genoFormat)
    
    alleles = genomeSite.alleles()
    
    if len(alleles) == 2:
        counts = genomeSite.asList(mode="count", samples=samples, missing = 9)
        
        genoOut.write("".join([str(c) for c in counts]) + "\n")
        
        if siteData["scaffold"] != scaf:
            #different scaffold from the last site
            #if using cumulative positions, change the offset for the last chrom
            if chrom is not None and args.cumulativePos: chromOffset[chrom] = pos
            #now get new scaf, chrom and pos
            scaf = siteData["scaffold"]
            try: chrom = chromDict[scaf]
            except: chrom = str(args.nullChrom)
        
        pos = siteData["position"] if not args.cumulativePos else siteData["position"] + chromOffset[chrom]
        
        snpOut.write("\t".join([str(linesDone), chrom, "0.0", str(pos), alleles[0], alleles[1]]) + "\n")
    
    linesDone += 1
    if linesDone % 100000 == 0: print linesDone, "lines done..."
    

genoOut.close
snpOut.close()
