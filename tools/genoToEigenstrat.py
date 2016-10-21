import sys,argparse,gzip
import genomics

chromInt = {"chr1":"1", "chr2":"2", "chr3":"3", "chr4":"4", "chr5":"5", "chr6":"6", "chr7":"7", "chr8":"8", "chr9":"9", "chr10":"10", "chr11":"11", "chr12":"12", "chr13":"13", "chr14":"14", "chr15":"15", "chr16":"16", "chr17":"17", "chr18":"18", "chr19":"19", "chr20":"20", "chrZ":"23", "chr21":"23"}

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Genotype format [otherwise will be inferred (slower)]", action = "store",choices = ["phased","diplo","paired"])
parser.add_argument("--genoOutFile", help="Output Eigenstrat geno file", action = "store", required = True)
parser.add_argument("--snpOutFile", help="Output Eigenstrat snp file", action = "store", required = True)
parser.add_argument("-s", "--samples", help="Analysis threads", action = "store")
parser.add_argument("--skipChecks", help="Skip genotype checks to speed things up", action = "store_true")

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

linesDone = 0

for line in genoFile:
    site = genomics.parseGenoLine(line)
    
    genomeSite = genomics.GenomeSite(genotypes = site.GTs, sampleNames = allNames, genoFormat=args.genoFormat, skipChecks = args.skipChecks)
    
    if len(genomeSite.alleles()) == 2:
        counts = genomeSite.asList(mode="count", samples=samples, missing = 9)

        genoOut.write("".join([str(c) for c in counts]) + "\n")
        
        snpOut.write(str(linesDone) + "\t" + chromInt[site.scaffold] + "\t" + "0.0" + "\t" + str(site.position) + "\n")

    linesDone += 1
    if linesDone % 100000 == 0: print linesDone, "lines done..."
    

genoOut.close
snpOut.close()
