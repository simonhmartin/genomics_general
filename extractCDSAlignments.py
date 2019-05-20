'''
script to extract CDS alignments given a gff file and a geno file.
It assumes the geno file is in "phased" format (eg T/T A|T). And it splits the sequences into two.
To output a single sequence for each geno file add argument --dontSplit
'''

import genomics
import sys, argparse, gzip, subprocess
from collections import defaultdict

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
parser.add_argument('--split', help="Split sequences into haplotypes", dest='split', action='store_true')
parser.add_argument('--no-split', help="Do not split sequences", dest='split', action='store_false')
#define ploidy if not 2
parser.add_argument("--ploidy", help="Ploidy for splitting phased sequences", action = "store", type=int, nargs="+", default=2)
parser.set_defaults(split=True)

args = parser.parse_args()

with gzip.open(args.gff,"r") if args.gff.endswith(".gz") else open(args.gff,"r") as gff: gffLines = gff.readlines()

sys.stderr.write("Parsing gene data from gff\n")
geneData = genomics.parseGenes(gffLines)

if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"w") if args.outFile.endswith(".gz") else open(args.outFile,"w")

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
        
        #for each exon we exytract the sequence for all individuals
        #and build a dictionary of position and genotype for each individual
        for i in range(geneData[scaffold][mRNA]["exons"]):
            
            #extract the sequence data for this exon
            
            genoStream = tabixStream(args.genoFile,
                                     chrom = scaffold,
                                     start = geneData[scaffold][mRNA]["cdsStarts"][i],
                                     end = geneData[scaffold][mRNA]["cdsEnds"][i], header=True)
            
            reader = genomics.GenoFileReader(genoStream, splitPhased=args.split, ploidy=args.ploidy)
            
            #if this is the first exon, initialise dictionary of genotype for each pos for each seqName
            if i == 0:
                nSeqs = len(reader.names)
                empty = ["N"]*nSeqs
                siteGTdict = defaultdict(lambda: empty)
            
            #now update if there is data for this window
            for siteData in reader.siteBySite(asDict=True if args.samples else False):
                GTs = [siteData["GTs"][name] for name in args.seqNames] if args.seqNames else siteData["GTs"]
                siteGTdict[siteData["position"]] = [genomics.complement(gt) for gt in GTs] if geneData[scaffold][mRNA]['strand'] == "-" else GTs
        
        cdsPositions = genomics.CDSpositions(geneData[scaffold][mRNA]['cdsStarts'],
                                             geneData[scaffold][mRNA]['cdsEnds'],
                                             geneData[scaffold][mRNA]['strand'])
        
        cdsSeqs = [[siteGTdict[position][i] for position in cdsPositions] for i in range(nSeqs)]
        
        if args.outFormat == "fasta":
            outputNames = [name + "_" + mRNA + " " + scaffold + " " +
                           str(geneData[scaffold][mRNA]['start']) + "-" + str(geneData[scaffold][mRNA]['end']) for name in reader.names]
        else: outputNames = [name + "_" + mRNA for name in reader.names]
        
        alnString = genomics.makeAlnString(outputNames,cdsSeqs,outFormat=args.outFormat, lineLen=None)
        outFile.write(alnString + "\n")

outFile.close()
