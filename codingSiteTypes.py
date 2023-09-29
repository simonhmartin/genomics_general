#!/usr/bin/env python

# script to annotate each coding site in the genome according to its type
# (codon position, synonymous / non-synonymous (variants only) and degeneracy.
# A vcf can be included to improve accuracy
# NOTE includion of a vcf requires that the vcf is tabix-indexed and tabix is installed. See https://github.com/samtools


import genomics
import sys, argparse, gzip, subprocess, itertools

################################################################################

def tabixStream(fileName, region = None, chrom = None, start=None, end=None, header=False):
    if not region: region = chrom+":"+str(start)+"-"+str(end)  
    
    if header:
        p = subprocess.Popen(['tabix', '-h', fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    else:
        p = subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True)
    
    return iter(p.communicate()[0].strip().split("\n"))

################################################################################

parser=argparse.ArgumentParser()

parser.add_argument("-a", "--annotation", help="Input annotation file (gff or gtf)", action = "store", required = True)
parser.add_argument("-f", "--format", help="Annotation file format", action = "store", choices = ("gff3", "gtf"), default="gff3")
parser.add_argument("-o", "--outFile", help="Output siteTypes file", action = "store")
parser.add_argument("-v", "--vcf", help="VCF file with variants (only ALTs will be considered)", action='store')
parser.add_argument("-r", "--ref", help="Genome .fa file. Must also have ", action='store', required = True)
parser.add_argument("--ignoreConflicts", help="Don't fail if two annotations give conflicting information about the same site", action='store_true')
parser.add_argument("--scaffoldLookup", help="Table of scaffold names (col1:reference and col2:annotation) if they differ", action = "store")
parser.add_argument("--useAnnotationScaffoldNames", help="Use the scaffold names in the annotation rather than the reference", action = "store_true")

args = parser.parse_args()

################################################################################

#get gene data
sys.stderr.write("Parsing annotation\n")
with gzip.open(args.annotation,"rt") if args.annotation.endswith(".gz") else open(args.annotation,"rt") as ann:
    geneData = genomics.parseGenes(ann.readlines(), fmt=args.format) 

#get scaffold names
sys.stderr.write("Loading reference genome\n")
with gzip.open(args.ref,"rt") if args.ref.endswith(".gz") else open(args.ref,"rt") as ref:
    scaffolds,_sequences_ = genomics.parseFasta(ref.read(), makeUppercase=True)
    sequences = {}
    for i,scaffold in enumerate(scaffolds): sequences[scaffold] = _sequences_[i]

#if a scaffold lookup table is provided, we need to change the names
if args.scaffoldLookup and args.useAnnotationScaffoldNames:
    #if we are going to use the annotation scaffold names, change the reference
    with open(args.scaffoldLookup) as lookup:
        scafNamesDict = dict([line.split() for line in lookup])
    sequences_renamed = {}
    scaffolds_renamed = []
    for scaffold in scaffolds:
        if scaffold in scafNamesDict.keys():
            sequences_renamed[scafNamesDict[scaffold]] = sequences[scaffold]
            scaffolds_renamed.append(scafNamesDict[scaffold])
        else:
            sys.stderr.write(f" WARNING!: {scaffold} is not in scaffoldLookup and will not be analysed\n")
    sequences = sequences_renamed
    scaffolds = scaffolds_renamed

if args.scaffoldLookup and args.useAnnotationScaffoldNames:
    #if we are going to use the reference scaffold names, change the gene data
    with open(args.scaffoldLookup) as lookup:
        scafNamesDict = dict([line.split()[::-1] for line in lookup])
    geneData_renamed = {}
    for scaffold in scaffolds: geneData_renamed[scaffold] = geneData[scafNamesDict[scaffold]]
    geneData = geneData_renamed

#open output
if not args.outFile: outFile = sys.stdout
else: outFile = gzip.open(args.outFile,"wt") if args.outFile.endswith(".gz") else open(args.outFile,"wt")

outFile.write("\t".join(["scaffold","position","codon_position","substitution_type","degeneracy"]) + "\n")

for scaffold in scaffolds:
    posData = {}
    
    positionsAnalysed = set()
    
    try:
        mRNAs = geneData[scaffold].keys()
        sys.stderr.write("Analysing {}: {} mRNAs\n".format(scaffold,len(mRNAs)))
    except:
        sys.stderr.write("Skipping {}. No annotated mRNAs\n".format(scaffold))
        continue
    
    counter = 0
    
    for mRNA in mRNAs:
        region = scaffold + ":" + str(geneData[scaffold][mRNA]["start"]) + "-" + str(geneData[scaffold][mRNA]["end"])
        sys.stderr.write("    Analysing mRNA {}: {}, {} exons\n".format(mRNA, region, geneData[scaffold][mRNA]["exons"]))
        
        siteAlleles = {}
        
        for i in range(geneData[scaffold][mRNA]["exons"]):
            
            start = geneData[scaffold][mRNA]["cdsStarts"][i]
            end = geneData[scaffold][mRNA]["cdsEnds"][i]
            
            #extract the sequence data for this gene, to be used below
            #first get sequence from reference
            
            siteAlleles.update(dict(zip(range(start,end+1), [set(base) for base in sequences[scaffold][start-1:end]])))
            
            #then, get variants
            if args.vcf:
                vcfStream = tabixStream(args.vcf,
                                        chrom = scaffold,
                                        start = geneData[scaffold][mRNA]["cdsStarts"][i],
                                        end = geneData[scaffold][mRNA]["cdsEnds"][i])
                
                for line in vcfStream:
                    if line != "":
                        CHROM,POS,ID,REF,ALT = line.split()[:5]
                        for a in ALT:
                            if a == "A" or a == "C"  or a == "G"  or a == "T" : siteAlleles[int(POS)].add(a)
        
        #remove any Ns
        for alleleSet in siteAlleles:
            try: alleleSet.remove("N")
            except: pass
        
        #for each CDS, extract a list of scaffold positions
        cdsPositions = genomics.CDSpositions(geneData[scaffold][mRNA]['cdsStarts'],
                                             geneData[scaffold][mRNA]['cdsEnds'],
                                             geneData[scaffold][mRNA]['strand'], trim=True)
        
        #make a list of lists for each codon. Each of which is a list of alleles
        codonAlleles = [[siteAlleles[cdsPositions[y]] if geneData[scaffold][mRNA]["strand"] == "+"
                         else genomics.complement(siteAlleles[cdsPositions[y]]) for y in range(x,x+3)] for x in range(len(cdsPositions))[::3]]
        
        #dictionary of substitution type and degeneracy for each position
        _posData_ = dict(zip(cdsPositions, [x for _codonAlleles_ in codonAlleles
                                            for x in zip(range(1,4),
                                                         genomics.synNon(*_codonAlleles_),
                                                         genomics.degeneracy(*_codonAlleles_))]))
        
        #now check whether we don't already have site type data for these seites from a previous gene
        previouslyAnalysedPositions = positionsAnalysed.intersection(cdsPositions)
        
        #check for conflicts if any sites are included in multiple RNAs
        for pos in previouslyAnalysedPositions:
            if posData[pos] != _posData_[pos]:
                if args.ignoreConflicts:
                    #print("WARNING: Position {} of {} occurs in two mRNAs giving conflicting site classifications. Setting output to NA\n".format(pos, scaffold), file=sys.stderr)
                    _posData_[pos] = ("NA","NA","NA")
                else:
                    raise AssertionError("Position {} of {} occurs in two mRNAs giving conflicting site classifications.\n".format(pos, scaffold))
        
        #update position data and the list of positions analysed 
        posData.update(_posData_)
        positionsAnalysed.update(cdsPositions)
        
        counter += 1
    
    # report number of genes for this scaffold
    sys.stderr.write("    Done analysing {} mRNAs. Writing output for {}\n".format(counter, scaffold))
    
    for pos in sorted(positionsAnalysed):
        outFile.write("\t".join([scaffold, str(pos), "\t".join([str(x) for x in posData[pos]])]) + "\n")
    
outFile.close()
