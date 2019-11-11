#!/usr/bin/env python

#Various functions to parse and filter genotype data from vcf files.
#If run independently it will pipe "genotypes" format to stdout. 

import argparse, sys, gzip, re, subprocess

import numpy as np

def GTtype(alleles):
    alleleSet = set(alleles)
    if len(alleleSet) > 1: return "Het"
    elif "0" in alleleSet: return "HomRef"
    elif "." in alleleSet: return "Missing"
    else: return "HomAlt"


class VcfSite:
    
    __slots__ = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_ALT', 'QUAL', 'FILTER', 'INFO', 'sampleNames', 'genoData', "alleleDict"]
    
    def __init__(self, elements=None, line=None, headers=None, headerLine=None, precompGenoData=None):
        assert((elements != None or line != None) and (headers != None or headerLine != None))
        if not headers: headers = headerLine.split()
        if not elements: elements = line.split()
        
        lineDict = dict(zip(headers,elements))
        
        self.CHROM = lineDict["#CHROM"]
        self.POS = int(lineDict["POS"])
        self.ID = lineDict["ID"]
        self.REF = lineDict["REF"]
        self.ALT = lineDict["ALT"].split(",")
        self.alleleDict = dict(zip([str(i) for i in range(len(self.ALT)+1)], [self.REF] + self.ALT))
        self.QUAL = lineDict["QUAL"]
        self.FILTER = lineDict["FILTER"]
        self.INFO = lineDict["INFO"].split(";")
        
        genoInfoNames = lineDict["FORMAT"].split(":")
        
        self.sampleNames = headers[9:]
        
        self.genoData = {}
        for sampleName in self.sampleNames:
            #if pre-compiled genotype data are available, try using those 
            if precompGenoData and lineDict[sampleName] in precompGenoData:
                self.genoData[sampleName] = precompGenoData[lineDict[sampleName]]
            else:
                #otherwise make dictionary for this sample
                self.genoData[sampleName] = dict(zip(genoInfoNames, lineDict[sampleName].split(":")))
                if "GT" in self.genoData[sampleName]:
                    self.genoData[sampleName]["alleles"] = tuple(self.genoData[sampleName]["GT"])[::2]
                    self.genoData[sampleName]["phase"] = "|" if "|" in self.genoData[sampleName]["GT"] else "/"
                if precompGenoData["__counter__"] < precompGenoData["__maxSize__"]:
                    precompGenoData[lineDict[sampleName]] = self.genoData[sampleName]
                    precompGenoData["__counter__"] += 1
    
    
    def getGenotype(self, sample, gtFilters = [], withPhase=True, asNumbers = False, missing = None, allowOnly=None, keepPartial=False, ploidy=None):
        genoData = self.genoData[sample]
        if missing is None:
            if asNumbers: missing = "."
            else: missing = "N"
        
        #check each gt filter
        passed = True
        for gtFilter in gtFilters:
            #first check that it's applicable
            if ("siteTypes" in gtFilter and self.getType() not in gtFilter["siteTypes"]): continue
            if ("gtTypes" in gtFilter and GTtype(genoData["alleles"]) not in gtFilter["gtTypes"]): continue
            if ("samples" in gtFilter and sample not in gtFilter["samples"]): continue
            #now check that it passes
            #might be a single value, nut could be several separated by commas. So will split in case
            values = np.array(genoData[gtFilter["flag"]].split(","), dtype=float)
            passed = np.all(gtFilter["min"] <= values) and np.all(values <= gtFilter["max"])
            #try: passed = gtFilter["min"] <= float(genoData[gtFilter["flag"]]) <= gtFilter["max"]
            #except: passed = False
            if not passed: break
        
        if ploidy is None: ploidy=len(genoData["alleles"])
        
        if passed:
            if not asNumbers:
                try:
                    sampleAlleles = [self.alleleDict[a] for a in genoData["alleles"]]
                    if allowOnly: sampleAlleles = [a if a in allowOnly else missing for a in sampleAlleles]
                    if not keepPartial: sampleAlleles = sampleAlleles if missing not in sampleAlleles else [missing]*ploidy
                
                except: sampleAlleles = [missing]*ploidy
            
            else:
                sampleAlleles = genoData["alleles"][:]
        
        
        else: sampleAlleles = [missing]*ploidy
        
        if withPhase: return genoData["phase"].join(sampleAlleles)
        else: return "".join(sampleAlleles)
    
    
    def getGenotypes(self, gtFilters = [], asList = False, withPhase=True, asNumbers = False,
                     samples = None, missing = None, allowOnly=None, keepPartial=False, ploidyDict=None):
        
        if not samples: samples = self.sampleNames
        output = {}
        for sample in samples:
            ploidy = ploidyDict[sample] if ploidyDict is not None else None
            output[sample] = self.getGenotype(sample, gtFilters=gtFilters, withPhase=withPhase, asNumbers=asNumbers,
                                              missing=missing, allowOnly=allowOnly, keepPartial=keepPartial, ploidy=ploidy)
        
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
            try: fields.append(self.genoData[sample][field])
            except: fields.append(missing)
        return fields


def parseHeaderLines(fileObj):
    output = {}
    output["contigs"] = []
    output["contigLengths"] = {}
    for line in fileObj:
        if line.startswith("##contig"):
            contigDataDict = dict([x.split("=") for x in re.split('<|>', line)[1].split(",")])
            elements = re.split('=|,|>', line)
            output["contigs"].append(contigDataDict["ID"])
            try: output["contigLengths"][contigDataDict["ID"]] = int(contigDataDict["length"])
            except: output["contigLengths"][contigDataDict["ID"]] = None
        
        if line.startswith("#CHROM"):
            output["mainHead"] = line
            elements = line.split()
            output["sampleNames"] = line.split()[9:]
            output["nSamples"] = len(output["sampleNames"])
            output["mainHeaders"] = elements
            break
    
    return output


def getHeadData(fileName):
    with gzip.open(fileName, "rt") if fileName.endswith(".gz") else open(fileName, "rt") as fileObj:
        return parseHeaderLines(fileObj)


def parseVcfSites(lines, mainHeaders, precomp=True, precompMaxSize=10000, excludeDuplicates=False):
    if precomp:
        precompGenoData = {}
        precompGenoData["__maxSize__"] = precompMaxSize
        precompGenoData["__counter__"] = 0
    else: precompGenoData = None
    
    if excludeDuplicates: lastChrom = lastPos = None
    
    for elements in lines:
        if isinstance(elements, str): elements = elements.split()
        if elements[0][0] == "#": continue
        if excludeDuplicates:
            if elements[0] == lastChrom and elements[1] == lastPos: continue
            lastChrom = elements[0]
            lastPos = elements[1]
        yield VcfSite(elements=elements, headers=mainHeaders, precompGenoData=precompGenoData)

def canFloat(string):
    try: float(string)
    except: return False
    return True

def parseGenotypeFilterArg(arg):
    try:
        gtfDict = dict([tuple(i.split("=")) for i in arg])
        for key in gtfDict.keys():
            assert key in ["flag","min","max", "siteTypes", "gtTypes", "samples"]
        for key in ["siteTypes", "gtTypes", "samples"]:
            if key in gtfDict: gtfDict[key] = gtfDict[key].split(",")
        gtfDict["min"] = float(gtfDict["min"]) if "min" in gtfDict else -np.inf
        gtfDict["max"] = float(gtfDict["max"]) if "max" in gtfDict else np.inf
        return gtfDict
    except: raise ValueError("Bad genotype filter specification. See help.")

###############################################################################################################
if __name__ == "__main__":


    ### parse arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inFile", help="Input vcf file", action = "store")
    parser.add_argument("-o", "--outFile", help="Output csv file", action = "store")

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
    
    parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+", default=[2])
    parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
    
    parser.add_argument("--field", help="Optional - format field to extract", action = "store")
    parser.add_argument("--missing", help="Value to use for missing data", action = "store")
    parser.add_argument("--outSep", help="Output separator", action = "store", default = "\t")

    args = parser.parse_args()

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
        sys.stderr.write("{} contigs will be included.".format(len(include)))
    
    if len(exclude) >= 1:
        exclude = set(exclude)
        sys.stderr.write("{} contigs will be excluded.".format(len(exclude)))
    
    gtFilters = [parseGenotypeFilterArg(gtf) for gtf in args.gtf] if args.gtf else []
    
    ##########################################################################################################################

    ### open files

    if args.inFile: inFile = gzip.open(args.inFile, "rt") if args.inFile.endswith(".gz") else open(args.inFile, "rt")
    else: inFile = sys.stdin


    if args.outFile: outFile = gzip.open(args.outFile, "w") if args.outFile.endswith(".gz") else open(args.outFile, "w")
    else: outFile = sys.stdout
    
    #header data
    headData = parseHeaderLines(inFile)
    
    #check specified samples are in first file. Otherwise use this entire set    
    if samples:
        for sample in samples: assert sample in headData["sampleNames"], "Specified sample name not in VCF header."
    else: samples = headData["sampleNames"]
    
    if args.ploidyFile is not None:
        with open(args.ploidyFile, "rt") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
    else:
        ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(samples)
        assert len(ploidy) == len(samples), "Incorrect number of ploidy values supplied."
        ploidyDict = dict(zip(samples,ploidy))


    ##########################################################################################################################

    outFile.write(args.outSep.join(["#CHROM", "POS"] + samples) + "\n")
    
    for vcfSite in parseVcfSites(inFile, headData["mainHeaders"]):
        if (exclude and vcfSite.CHROM in exclude) or (include and vcfSite.CHROM not in include): continue
        if args.skipMono and vcfSite.getType() is "mono": continue
        if args.minQual and canFloat(vcfSite.QUAL) and float(vcfSite.QUAL) < args.minQual: continue
        if args.field is not None: output = vcfSite.getGenoField(args.field,samples=samples, missing=args.missing)
        else:
            allowed=["A","C","G","T"] if args.skipIndels else None
            output = vcfSite.getGenotypes(gtFilters,asList=True,withPhase=True,samples=samples,missing=args.missing,
                                            allowOnly=allowed,keepPartial=False,ploidyDict=ploidyDict)
        outFile.write(args.outSep.join([vcfSite.CHROM, str(vcfSite.POS)] + output) + "\n")
    
    outFile.close()
