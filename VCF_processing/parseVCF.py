#!/usr/bin/env python

#Various functions to parse and filter genotype data from vcf files.
#If run independently it will pipe "genotypes" format to stdout. 

import argparse, sys, gzip, re, subprocess


class vcfGenoData:
    
    def __init__(self, formatList, genoList):
        genoInfoNames = formatList
        genoInfo = genoList
        for x in range(len(genoInfo)):
            setattr(self, genoInfoNames[x], genoInfo[x])
        if hasattr(self,"GT"):
            if len(self.GT) == 3:
                self.phase = self.GT[1]
                self.alleles = self.GT.split(self.phase)
            elif len(self.GT) == 1:
                self.phase = ""
                self.alleles = self.GT
            else: raise ValueError, "Error parsing genotype. Check genotype field."
    
    def getType(self):
        if not hasattr(self,"GT"): return None
        elif len(self.alleles) == 1: return "Haploid"
        elif self.alleles[0] == self.alleles[1] == "0": return "HomRef"
        elif self.alleles[0] != self.alleles[1]: return "Het"
        elif self.alleles[0] == self.alleles[1] != ".": return "HomAlt"
        else: return "Missing"


class VcfSite:
    
    def __init__(self, Line, head):
        headers = head.split()
        elements = Line.split()
        assert len(headers) == len(elements), "Error - vcf line has different number of elements to header"
        
        lineDict = dict(zip(headers,elements))
        
        self.CHROM = lineDict["#CHROM"]
        self.POS = int(lineDict["POS"])
        self.ID = lineDict["ID"]
        self.REF = lineDict["REF"]
        self.ALT = lineDict["ALT"].split(",")
        self.QUAL = lineDict["QUAL"]
        self.FILTER = lineDict["FILTER"]
        self.INFO = lineDict["INFO"].split(";")
                
        genoInfoNames = lineDict["FORMAT"].split(":")
        
        self.sampleNames = headers[9:]
        
        self.genoData = {}
        for sampleName in self.sampleNames:
            self.genoData[sampleName] = vcfGenoData(genoInfoNames, lineDict[sampleName].split(":"))
    
    
    def getGenotype(self, sample, gtFilters = [], withPhase=True, asNumbers = False, missing = None, allowOnly=None):
        genoData = self.genoData[sample]
        if missing is None:
            if asNumbers: missing = "."
            else: missing = "N"
        
        #check each gt filter
        passed = True
        for gtFilter in gtFilters:
            #first check that it's applicable
            if (gtFilter.siteTypes and self.getType() not in gtFilter.siteTypes) or (gtFilter.gtTypes and genoData.getType() not in gtFilter.gtTypes) or (gtFilter.samples and sample not in gtFilter.samples):
                continue
            #now check that it passes
            try:
                assert gtFilter.Min <= float(getattr(genoData, gtFilter.flag)) <= gtFilter.Max
            except:
                passed = False
                break
        
        ploidy=len(genoData.alleles)

        if passed:
            sampleAlleles = genoData.alleles
            if not asNumbers:
                alleles = [self.REF]+self.ALT
                if allowOnly: alleles = [a if a in allowOnly else missing for a in alleles] 
                try: sampleAlleles = [alleles[int(a)] for a in sampleAlleles]
                except: sampleAlleles = [missing]*ploidy
        
        else: sampleAlleles = [missing]*ploidy
        
        if withPhase: return "/".join(sampleAlleles)
        else: return "".join(sampleAlleles)

    
    def getGenotypes(self, gtFilters = [], asList = False, withPhase=True, asNumbers = False, samples = None, missing = None, allowOnly=None):
        if not samples: samples = self.sampleNames
        output = {}
        for sample in samples:
            output[sample] = self.getGenotype(sample, gtFilters=gtFilters, withPhase=withPhase, asNumbers=asNumbers, missing=missing, allowOnly=allowOnly)
        
        if asList: return [output[sample] for sample in samples]
        
        return output
    
    
    def getType(self):
        if len(self.REF) == 1:
            if self.ALT == ["."]: return "mono"
            elif max([len(a) for a in self.ALT]) == 1: return "SNP"
            else: return "indel"
        else: return "indel"
    
    def getGenoField(self, field, samples = None, missing=None):
        if missing is None: missing = "."
        if samples is None: samples = self.sampleNames
        fields = []
        for sample in samples:
            try: fields.append(getattr(self.genoData[sample], field))
            except: fields.append(missing)
        return fields

class gtFilter:
    
    def __init__(self, flag, Min = None, Max = None, samples = None, siteTypes = None, gtTypes = None):
        if not Min: Min = "-inf"
        Min = float(Min)
        if not Max: Max = "inf"
        Max = float(Max)
        self.flag = flag
        self.Min = Min
        self.Max = Max
        self.samples = samples
        self.siteTypes = siteTypes
        self.gtTypes = gtTypes     



class HeadData:
    
    def __init__(self, headerLines):
        self.contigs = []
        self.contigLengths = {}
        for Line in headerLines:
            if Line[:8] == "##contig":
                try:
                    elements = re.split('=|,|>', Line)
                    self.contigs.append(elements[2])
                    try: self.contigLengths[elements[2]] = int(elements[4])
                    except: self.contigLengths[elements[2]] = None
                    
                except:
                    pass
            
            if Line[:6] == "#CHROM":
                self.sampleNames = Line.split()[9:]
                self.mainHead = Line
            

def getHeadData(fileName):
    headLines = []
    with gzip.open(fileName, "r") if fileName.endswith(".gz") else open(fileName, "r") as fileObj:
        for line in fileObj:
            if line.startswith("#"): headLines.append(line)
            else: break
    return HeadData(headLines)

class Reader:
    
    def __init__(self, fileObj):
        self.fileObj = fileObj
        headLines = []
        line = fileObj.readline()
        while(line[:1] == "#"):
            headLines.append(line)
            if line[:6] != "#CHROM":
                line = fileObj.readline()
            else:
                break
        headData = HeadData(headLines)
        self.contigs = headData.contigs
        self.contigLengths = headData.contigLengths
        self.sampleNames = headData.sampleNames
        self.mainHead = headData.mainHead
    
    def lines(self):
        line = self.fileObj.readline()
        while len(line) >= 1:
            yield line
            line = self.fileObj.readline()
    
    def sites(self):
        line = self.fileObj.readline()
        while len(line) >= 1:
            site = VcfSite(line, self.mainHead)
            yield site
            line = self.fileObj.readline()

def tabixStream(fileName, chrom, start, end):
    region = chrom+":"+str(start)+"-"+str(end)  
    return subprocess.Popen(['tabix',fileName, region], stdout=subprocess.PIPE, bufsize=1)

def tabixLines(fileName, chrom, start, end):
    stream = tabixStream(fileName, chrom, start, end)
    for line in iter(stream.stdout.readline, ""): yield line

def tabixSites(fileName, chrom, start, end, mainHead):
    lineGen = tabixLines(fileName, chrom, start, end)
    for line in lineGen: yield VcfSite(line, mainHead)

def canFloat(string):
    try: float(string)
    except: return False
    return True


###############################################################################################################
if __name__ == "__main__":


    ### parse arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input vcf file", action = "store")
    parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")

    #specific samples
    parser.add_argument("-s", "--samples", help="sample names (separated by commas)", action='store')

    #contigs
    parser.add_argument("--include", help="include contigs (separated by commas)", action='store')
    parser.add_argument("--includeFile", help="File of contigs (one per line)", action='store')
    parser.add_argument("--exclude", help="exclude contigs (separated by commas)", action='store')
    parser.add_argument("--excludeFile", help="File of contigs (one per line)", action='store')


    #vcf parsing arguments
    parser.add_argument("--minQual", help="Minimum QUAL for a site", type=int, action = "store")
    parser.add_argument("--gtf", help="Genotype filter. Syntax: flag=X min=X max=X siteTypes=X,X.. gtTypes=X,X.. samples=X,X..", action = "append", nargs = '+')

    parser.add_argument("--skipIndels", help="Skip indels", action = "store_true")
    parser.add_argument("--skipMono", help="Skip monomorphic sites", action = "store_true")
    
    parser.add_argument("--field", help="Optional - format field to extract", action = "store")
    parser.add_argument("--missing", help="Value to use for missing data", action = "store")
    parser.add_argument("--outSep", help="Output separator", action = "store", default = "\t")

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    samples = args.samples

    if samples: samples = samples.split(",")

    include = []
    exclude = []

    if args.include: include += args.include.split(",")
    if args.exclude: exclude += args.exclude.split(",")

    if args.includeFile:
        with open(args.includeFile, 'r') as includeFile:
            include += [c.strip() for c in includeFile.read().split("\n")]

    if args.excludeFile:
        with open(args.excludeFile, 'r') as excludeFile:
            exclude += [c.strip() for c in excludeFile.read().split("\n")]

    if len(include) >= 1:
        include = set(include)
        print >> sys.stderr, len(include), "contigs will be included."
        
    if len(exclude) >= 1:
        exclude = set(exclude)
        print >> sys.stderr, len(exclude), "contigs will be excluded."

    minQual = args.minQual

    gtFilters = []
    if args.gtf:
        for gtf in args.gtf:
            try:
                gtfDict = dict([tuple(i.split("=")) for i in gtf])
                for key in gtfDict.keys():
                    assert key in ["flag","min","max", "siteTypes", "gtTypes", "samples"]
                for x in ["siteTypes", "gtTypes", "samples"]:
                    if x in gtfDict.keys(): gtfDict[x] = gtfDict[x].split(",")
                    else: gtfDict[x] = None
                for x in ["min", "max"]:
                    if x not in gtfDict.keys(): gtfDict[x] = None
                gtFilters.append(gtFilter(gtfDict["flag"], gtfDict["min"], gtfDict["max"], gtfDict["samples"], gtfDict["siteTypes"], gtfDict["gtTypes"]))
            except:
                print >> sys.stderr, "Bad genotype filter specification. See help."  
                raise


    skipIndels = args.skipIndels
    skipMono = args.skipMono

    outSep = args.outSep
    ##########################################################################################################################

    ### open files

    if infile:
        if infile[-3:] == ".gz":
            In = gzip.open(infile, "r")
        else:
            In = open(infile, "r")
    else:
        In = sys.stdin

    if outfile:
        if outfile[-3:] == ".gz":
            Out = gzip.open(outfile, "w")
        else:
            Out = open(outfile, "w")
    else:
        Out = sys.stdout


    ### read through header for all input files

    vcf = Reader(In)

    #check specified samples are in first file. Otherwise use this entire set

    if samples:
        for sample in samples: assert sample in vcf.sampleNames, "Specified sample name not in VCF header."
    else: samples = vcf.sampleNames

    ##########################################################################################################################

    Out.write(" ".join(["#CHROM", "POS"] + samples) + "\n")

    vcfSites = vcf.sites()

    for vcfSite in vcfSites:
        if (exclude and vcfSite.CHROM in exclude) or (include and vcfSite.CHROM not in include): continue
        #print >> sys.stderr, vcfSite.CHROM, vcfSite.POS, vcfSite.REF, vcfSite.ALT, vcfSite.getType()
        if skipIndels and vcfSite.getType() is "indel": continue
        if skipMono and vcfSite.getType() is "mono": continue
        if minQual and canFloat(vcfSite.QUAL) and float(vcfSite.QUAL) < minQual: continue
        if args.field is not None: output = vcfSite.getGenoField(args.field,samples=samples, missing=args.missing)
        else: output = vcfSite.getGenotypes(gtFilters,asList=True,withPhase=True,samples=samples,missing=args.missing,allowOnly="ACGT")
        Out.write(outSep.join([vcfSite.CHROM, str(vcfSite.POS)] + output) + "\n")

