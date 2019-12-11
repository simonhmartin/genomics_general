'''
script to extract CDS alignments given a gff file and a geno file.
It assumes the geno file is in "phased" format (eg T/T A|T). And it splits the sequences into two.
To output a single sequence for each geno file add argument --dontSplit
'''

import genomics
import sys, argparse, gzip, subprocess
from collections import defaultdict
import numpy as np

def tabixStream(fileName, region = None, chrom = None, start=None, end=None, header=False):
    if not region: region = chrom+":"+str(start)+"-"+str(end)  
    
    if header:
        p = subprocess.Popen(['tabix', '-h', fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    else:
        p = subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    
    return iter(p.communicate()[0].strip().split("\n"))

##################################################

parser=argparse.ArgumentParser()

parser.add_argument("--gff", help="Input gff file", action = "store", required = True)
parser.add_argument("-o", "--outFile", help="Output sequence file", action = "store")
parser.add_argument("--outFormat", help="Output file format", action = "store", choices = ["fasta","phylip"], default="phylip")
parser.add_argument("-g", "--genoFile", help="Genotype file", action='store', required = True)
parser.add_argument("-s", "--samples", help="Samples to include", nargs="+", action='store')
parser.add_argument("-t", "--targets", help="Target mRNA names", nargs="+", action='store')
parser.add_argument('--split', help="Split sequences into haplotypes", dest='split', action='store_true')
parser.add_argument('--no-split', help="Do not split sequences", dest='split', action='store_false')
#define ploidy if not 2
parser.add_argument("--ploidy", help="Ploidy for splitting phased sequences", action = "store", type=int, nargs="+", default=2)
parser.set_defaults(split=True)

args = parser.parse_args()

with gzip.open(args.gff,"rt") if args.gff.endswith(".gz") else open(args.gff,"rt") as gff: gffLines = gff.readlines()

sys.stderr.write("Parsing gene data from gff\n")
geneData = genomics.parseGenes(gffLines, targets=args.targets)

if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"wt") if args.outFile.endswith(".gz") else open(args.outFile,"wt")

###################################################

#extract each scaffold from the geno file, and the genes for each scaffold and write them out

for scaffold in geneData.keys():
    mRNAs = geneData[scaffold].keys()
    sys.stderr.write("Extracting " + str(len(mRNAs)) + " gene sequences from " + scaffold + "\n")
    for mRNA in mRNAs:
        region = scaffold + ":" + str(geneData[scaffold][mRNA]["start"]) + "-" + str(geneData[scaffold][mRNA]["end"])
        
        if geneData[scaffold][mRNA]["exons"] < 1:
            sys.stderr.write("    Skipping mRNA {}: {}. No exons\n".format(mRNA, region))
            continue
        
        sys.stderr.write("    Extracting mRNA {}: {}, {} exons\n".format(mRNA, region, geneData[scaffold][mRNA]["exons"]))
        
        #get the start and end positions of CDSs in order
        strand = geneData[scaffold][mRNA]['strand']
        order = np.argsort(geneData[scaffold][mRNA]["cdsStarts"])
        if strand == "-": order = order[::-1]
        cdsStarts = [geneData[scaffold][mRNA]["cdsStarts"][x] for x in order]
        cdsEnds = [geneData[scaffold][mRNA]["cdsEnds"][x] for x in order]
        
        #for each exon we exytract the sequence for all individuals
        #and build a dictionary of position and genotype for each individual
        for i in range(geneData[scaffold][mRNA]["exons"]):
            
            #extract the sequence data for this exon
            
            genoStream = tabixStream(args.genoFile,
                                     chrom = scaffold,
                                     start = cdsStarts[i],
                                     end = cdsEnds[i], header=True)
            
            reader = genomics.GenoFileReader(genoStream, splitPhased=args.split, ploidy=args.ploidy)
            
            #if this is the first exon, initialise dictionary of genotype for each pos for each seqName
            if i == 0:
                nSeqs = len(reader.names)
                empty = ["N"]*nSeqs
                siteGTdict = defaultdict(lambda: empty)
            
            #now update if there is data for this window
            for siteData in reader.siteBySite(asDict=True if args.samples else False):
                GTs = [siteData["GTs"][name] for name in args.samples] if args.samples else siteData["GTs"]
                siteGTdict[siteData["position"]] = [genomics.complement(gt) for gt in GTs] if strand == "-" else GTs
        
        cdsPositions = genomics.CDSpositions(cdsStarts,
                                             cdsEnds,
                                             strand)
        
        cdsSeqs = [[siteGTdict[position][i] for position in cdsPositions] for i in range(nSeqs)]
        
        if args.outFormat == "fasta":
            outputNames = [name + "_" + mRNA + " " + scaffold + " " +
                           str(geneData[scaffold][mRNA]['start']) + "-" + str(geneData[scaffold][mRNA]['end']) for name in reader.names]
        else: outputNames = [name + "_" + mRNA for name in reader.names]
        
        alnString = genomics.makeAlnString(outputNames,cdsSeqs,outFormat=args.outFormat, lineLen=None)
        outFile.write(alnString + "\n")

outFile.close()
