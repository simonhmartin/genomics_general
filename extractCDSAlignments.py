'''
script to extract CDS alignments given a gff file and a geno file.
It assumes the geno file is in "phased" format (eg T/T A|T). And it splits the sequences into two.
To output a single sequence for each geno file add argument --dontSplit
'''

import genomics
import sys, argparse, gzip, subprocess

##################################################

parser=argparse.ArgumentParser()

parser.add_argument("--gff", help="Input gff file", action = "store", required = True)
parser.add_argument("-o", "--outFile", help="Output gff file", action = "store")
parser.add_argument("-g", "--genoFile", help="Genotype file", action='store', required = True)
parser.add_argument("-s", "--samples", help="Samples to include", nargs="+", action='store')
parser.add_argument('--split', help="Split diploid sequences into two haplotypes", dest='split', action='store_true')
parser.add_argument('--no-split', help="Do not split sequences", dest='split', action='store_false')
parser.set_defaults(split=True)

args = parser.parse_args()
#args = parser.parse_args(["--gff", "/scratch/shm45/Hmel2/Hmel2.cortex.gff", "-g", "/zoo/disk1/shm45/vcf/rosina/ros10.Hmel2.bwa.default.HC.DP8.Hmel215006.geno.gz"])

with gzip.open(args.gff,"r") if args.gff.endswith(".gz") else open(args.gff,"r") as gff: gffLines = gff.readlines()

geneData = genomics.parseGenes(gffLines)

if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"w") if args.outFile.endswith(".gz") else open(args.outFile,"w")

###################################################

#extract each scaffold from the geno file, and the genes for each scaffold and write them out

for scaffold in geneData.keys():
    mRNAs = geneData[scaffold].keys()
    sys.stderr.write("Extracting " + str(len(mRNAs)) + " gene sequences from " + scaffold + "\n")
    for mRNA in mRNAs:
        sys.stderr.write(mRNA + "\n")
        region = scaffold + ":" + str(geneData[scaffold][mRNA]["start"]) + "-" + str(geneData[scaffold][mRNA]["end"])
        sys.stderr.write("Getting region " + region + " from geno file...\n") 
        genoStream = subprocess.Popen(['tabix', '-h', args.genoFile, region], stdout=subprocess.PIPE)
        window = genomics.parseGenoFile(genoStream.stdout, names=args.samples, includePositions=True, splitPhased = args.split)
        seqDict = window.seqDict()
        seqNames=seqDict.keys()
        sys.stderr.write("Extracting CDS...\n") 
        CDSseqs = [genomics.CDS(seqDict[name], window.positions, geneData[scaffold][mRNA]['cdsStarts'], geneData[scaffold][mRNA]['cdsEnds'],geneData[scaffold][mRNA]['strand']) for name in seqNames]
        
        outputNames = [name + "_" + mRNA + " " + scaffold + " " +
                       str(geneData[scaffold][mRNA]['start']) + "-" + str(geneData[scaffold][mRNA]['end']) for name in seqNames]
        fastaString = genomics.makeAlnString(outputNames,CDSseqs,format="fasta", lineLen=None)
        outFile.write(fastaString + "\n")

outFile.close()
