import sys,argparse,gzip
import genomics

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","paired"])
parser.add_argument("--genoOutFile", help="Output Eigenstrat geno file", action = "store", required = True)
parser.add_argument("--snpOutFile", help="Output Eigenstrat snp file", action = "store", required = True)
parser.add_argument("-s", "--samples", help="Analysis threads", action = "store")
parser.add_argument("--skipChecks", help="Skip genotype checks to speed things up", action = "store_true")
parser.add_argument("--chromFile", help="Input chromosome number file", action = "store")


args = parser.parse_args()

if args.genoFile:
    genoFile = gzip.open(args.genoFile,"r") if args.genoFile.endswith(".gz") else open(args.genoFile, "r")
else: genoFile = sys.stdin

genoOut = open(args.genoOutFile, "w")
snpOut = open(args.snpOutFile, "w")

headers = genoFile.readline().split()
allNames = headers[2:]

if args.samples is None: samples = allNames
else: samples = args.samples


if args.chromFile:
    chromsProvided = True
    chromDict = {}
    with open(args.chromFile, "r") as chromFile:
        for line in chromFile:
            scaf, chrom = line.split()
            chromDict[scaf] = chrom

else: chromsProvided = False

linesDone = 0

for line in genoFile:
    site = genomics.parseGenoLine(line)
    
    genomeSite = genomics.GenomeSite(genotypes = site.GTs, sampleNames = allNames, genoFormat=args.genoFormat, skipChecks = args.skipChecks)
    
    if len(genomeSite.alleles()) == 2:
        counts = genomeSite.asList(mode="count", samples=samples, missing = 9)

        genoOut.write("".join([str(c) for c in counts]) + "\n")
        
        if chromsProvided:
            try: chromNo = chromDict[site.scaffold]
            except: chromNo = "0"
            snpOut.write(str(linesDone) + "\t" + chromNo + "\t" + "0.0" + "\t" + str(site.position) + "\n")
        else:
            snpOut.write(str(linesDone) + "\t0\t" + "0.0" + "\t" + str(site.position) + "\n")

    linesDone += 1
    if linesDone % 100000 == 0: print linesDone, "lines done..."
    

genoOut.close
snpOut.close()
