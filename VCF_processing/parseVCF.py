#!/usr/bin/env python

#Various functions to parse and filter genotype data from vcf files.
#If run independently it will pipe "genotypes" format to stdout. 

import argparse, sys, gzip, re, subprocess

from collections import defaultdict

import numpy as np


def GTtype(alleles):
    alleleSet = set(alleles)
    if len(alleleSet) > 1: return "Het"
    elif "0" in alleleSet: return "HomRef"
    elif "." in alleleSet: return "Missing"
    else: return "HomAlt"


re_cigar = re.compile("\d+|[MXDI]")

re_phaser = re.compile("[/|]")

def simplifyAlt(alt, cigar, missing="N"):
    l = re_cigar.findall(cigar)
    i = 0
    simp = ""
    try:
        for x in range(0,len(l),2):
            label = l[x+1]
            n = int(l[x])
            if label == "M" or label == "X":
                #we have a matching stretch
                simp += alt[i:i + n]
                i+=n
            elif label == "I":
                #insertion, so just advance index
                i += n
            elif label == "D":
                #deletion, so add missing without advancing
                simp += "".join([missing]*n)
    except:
        raise ValueError("Malformed CIGAR: " + cigar)
    
    return(simp)


class VcfSite:
    
    __slots__ = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'REFlen', 'nALT', 'lenMatchDict', 'QUAL', 'FILTER', 'INFO', 'sampleNames', 'genoData', "alleleDict"]
    
    def __init__(self, elements=None, line=None, headers=None, headerLine=None, precompGenoData=None,
                 parseINFO=False, simplifyALT=False):
        assert((elements != None or line != None) and (headers != None or headerLine != None))
        if not headers: headers = headerLine.split()
        if not elements: elements = line.split()
        
        lineDict = dict(zip(headers,elements))
        
        self.CHROM = lineDict["#CHROM"]
        self.POS = int(lineDict["POS"])
        self.ID = lineDict["ID"]
        self.REF = lineDict["REF"]
        self.REFlen = len(self.REF)
        self.ALT = lineDict["ALT"].split(",") if lineDict["ALT"] != "." else []
        self.nALT = len(self.ALT)
        self.QUAL = lineDict["QUAL"]
        self.FILTER = lineDict["FILTER"]
        
        if parseINFO or simplifyALT:
            self.INFO = dict([x.split("=") for x in lineDict["INFO"].split(";")])
        
        if simplifyALT:
            cigars = self.INFO["CIGAR"].split(",")
            for x in range(self.nALT):
                self.ALT[x] = simplifyAlt(self.ALT[x], cigars[x])
        
        self.alleleDict = dict(zip([str(i) for i in range(self.nALT+1)], [self.REF] + self.ALT))
        self.lenMatchDict = dict([(key, len(self.alleleDict[key]) == self.REFlen) for key in self.alleleDict])
        
        genoInfoNames = lineDict["FORMAT"].split(":")
        
        self.sampleNames = headers[9:]
        
        self.genoData = {}
        for sampleName in self.sampleNames:
            #if pre-compiled genotype data are available, try using those 
            #need to check is both the format headings AND individual data match. So use a tuple as index
            if precompGenoData and (lineDict["FORMAT"], lineDict[sampleName],) in precompGenoData:
                self.genoData[sampleName] = precompGenoData[(lineDict["FORMAT"], lineDict[sampleName],)]
            else:
                #otherwise make dictionary for this sample
                self.genoData[sampleName] = dict(zip(genoInfoNames, lineDict[sampleName].split(":")))
                if "GT" in self.genoData[sampleName]:
                    self.genoData[sampleName]["alleles"] = tuple(re_phaser.split(self.genoData[sampleName]["GT"]))
                    self.genoData[sampleName]["phase"] = "|" if "|" in self.genoData[sampleName]["GT"] else "/"
                if precompGenoData["__counter__"] < precompGenoData["__maxSize__"]:
                    precompGenoData[(lineDict["FORMAT"], lineDict[sampleName],)] = self.genoData[sampleName]
                    precompGenoData["__counter__"] += 1
    
    def getSiteType(self):
        if len(self.ALT) == 0: return 'MONO'
        if np.all(list(self.lenMatchDict.values())): return 'SNP'
        return 'INDEL'
    
    def getGenotype(self, sample, gtFilters = [], withPhase=True, asNumbers = False, missing = None,
                    allowOnly=None, mustMatchREFlen=False, keepPartial=False, ploidy=None, ploidyMismatchToMissing=False, expandMulti=False):
        
        genoData = self.genoData[sample]
        if missing is None:
            if asNumbers: missing = "."
            else:
                missing = "N" if not expandMulti or self.REFlen==1 else ["N"]*self.REFlen
        
        #check each gt filter
        passed = True
        for gtFilter in gtFilters:
            #first check that it's applicable
            if ("siteTypes" in gtFilter and self.getSiteType() not in gtFilter["siteTypes"]): continue
            if ("gtTypes" in gtFilter and GTtype(genoData["alleles"]) not in gtFilter["gtTypes"]): continue
            if ("samples" in gtFilter and sample not in gtFilter["samples"]): continue
            #now check that it passes
            #might be a single value, nut could be several separated by commas. So will split in case
            try:
                values = np.array(genoData[gtFilter["flag"]].split(","), dtype=float)
                passed = np.all(gtFilter["min"] <= values) and np.all(values <= gtFilter["max"])
            except: passed = False
            #try: passed = gtFilter["min"] <= float(genoData[gtFilter["flag"]]) <= gtFilter["max"]
            #except: passed = False
            if not passed: break
        
        if ploidy is None: ploidy=len(genoData["alleles"])
        elif ploidy != len(genoData["alleles"]):
            if ploidyMismatchToMissing:
                passed=False
            else:
                raise ValueError("Sample {} at {}:{} genotype {} does not match explected ploidy of {}".format(sample, self.CHROM, self.POS,
                                                                                                           self.genoData[sample]["GT"], ploidy)) 
        
        if passed:
            if not asNumbers:
                try:
                    #retrieve alleles, but check if lengths must match and if they don't add a missing allele
                    sampleAlleles = [self.alleleDict[a] if (not mustMatchREFlen or self.lenMatchDict[a]) else missing for a in genoData["alleles"]]
                    if allowOnly: sampleAlleles = [a if a in allowOnly else missing for a in sampleAlleles]
                    if not keepPartial: sampleAlleles = sampleAlleles if missing not in sampleAlleles else [missing]*ploidy
                
                except:
                    sampleAlleles = [missing]*ploidy
            
            else:
                sampleAlleles = genoData["alleles"][:]
        
        
        else: sampleAlleles = [missing]*ploidy
        
        sep = genoData["phase"] if withPhase else ""
        
        if expandMulti:
            return tuple(sep.join([a[i] for a in sampleAlleles]) for i in range(self.REFlen))
        
        return sep.join(sampleAlleles)
    
    
    def getGenotypes(self, gtFilters = [], asList = False, withPhase=True, asNumbers = False,
                     samples = None, missing = None, allowOnly=None, mustMatchREFlen=False,
                     keepPartial=False, ploidyDict=None, ploidyMismatchToMissing=False, expandMulti=False):
        
        if not samples: samples = self.sampleNames
        output = {}
        for sample in samples:
            ploidy = ploidyDict[sample] if ploidyDict is not None else None
            output[sample] = self.getGenotype(sample, gtFilters=gtFilters, withPhase=withPhase, asNumbers=asNumbers,
                                              missing=missing, allowOnly=allowOnly, mustMatchREFlen=mustMatchREFlen,
                                              keepPartial=keepPartial, ploidy=ploidy, ploidyMismatchToMissing=ploidyMismatchToMissing,
                                              expandMulti=expandMulti)
        
        if asList: return [output[sample] for sample in samples]
        
        return output
    
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
            contigDataDict = dict([x.split("=", maxsplit=1) for x in re.split('<|>', line)[1].split(",")])
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


def parseVcfSites(lines, mainHeaders, precomp=True, precompMaxSize=10000, excludeDuplicates=False, parseINFO=False, simplifyALT=False):
    if precomp:
        precompGenoData = {}
        precompGenoData["__maxSize__"] = precompMaxSize
        precompGenoData["__counter__"] = 0
    else: precompGenoData = None
    
    if excludeDuplicates: lastChrom = lastPos = None
    
    for elements in lines:
        if isinstance(elements, str): elements = elements.split()
        if len(elements) == 0 or elements[0][0] == "#": continue
        if excludeDuplicates:
            if elements[0] == lastChrom and elements[1] == lastPos: continue
            lastChrom = elements[0]
            lastPos = elements[1]
        yield VcfSite(elements=elements, headers=mainHeaders, precompGenoData=precompGenoData, parseINFO=parseINFO, simplifyALT=simplifyALT)

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


def addArgs(parser, requireInfile=False):
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
    parser.add_argument("--excludeDuplicates", help="Only include the first in a series of duplicated positions", action = "store_true")
    parser.add_argument("--simplifyALT", help="Simplify multi-site alternate alleles using CIGAR (as in Freebayes output)", action = "store_true")
    parser.add_argument("--expandMulti", help="Expand multi-site alleles (also Sets simplifyALT to True)", action = "store_true")
    parser.add_argument("--maxREFlen", help="Maximum length for refernece allele", action = "store", type=int)
    
    parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, default=2)
    parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
    parser.add_argument("--ploidyMismatchToMissing", help="Set genotypes with mismatched ploidy to missing", action = "store_true")
    parser.add_argument("--keepPartial", help="Keep genotypes where some but not all alleles are missing", action = "store_true")
    parser.add_argument("--addRefTrack", help="Add a third column with the header REF and the reference allele", action = "store_true")
    
    parser.add_argument("--field", help="Optional - format field to extract", action = "store")
    parser.add_argument("--missing", help="Value to use for missing data", action = "store")
    parser.add_argument("--outSep", help="Output separator", action = "store", default = "\t")


def parseIncludeExcludeArgs(args):
    include = []
    exclude = []

    if args.include: include += args.include.split(",")
    if args.exclude: exclude += args.exclude.split(",")

    if args.includeFile:
        with open(args.includeFile, 'rt') as includeFile:
            include += [c.strip() for c in includeFile.read().split("\n")]

    if args.excludeFile:
        with open(args.excludeFile, 'rt') as excludeFile:
            exclude += [c.strip() for c in excludeFile.read().split("\n")]

    if len(include) >= 1:
        include = set(include)
        sys.stderr.write("{} contigs will be included.".format(len(include)))
    
    if len(exclude) >= 1:
        exclude = set(exclude)
        sys.stderr.write("{} contigs will be excluded.".format(len(exclude)))
    
    return (include,exclude,)


###############################################################################################################
if __name__ == "__main__":


    ### parse arguments

    parser = argparse.ArgumentParser()
    
    addArgs(parser)
    
    parser.add_argument("-i", "--inFile", help="Input vcf file", action = "store")
    
    args = parser.parse_args()
    
    samples = args.samples.split(",") if args.samples else None
    
    include,exclude = parseIncludeExcludeArgs(args)
    
    gtFilters = [parseGenotypeFilterArg(gtf) for gtf in args.gtf] if args.gtf else []
    
    simplifyALT = args.simplifyALT or args.expandMulti
    
    ##########################################################################################################################

    ### open files

    if args.inFile: inFile = gzip.open(args.inFile, "rt") if args.inFile.endswith(".gz") else open(args.inFile, "rt")
    else: inFile = sys.stdin


    if args.outFile: outFile = gzip.open(args.outFile, "wt") if args.outFile.endswith(".gz") else open(args.outFile, "wt")
    else: outFile = sys.stdout
    
    #header data
    headData = parseHeaderLines(inFile)
    
    #check specified samples are in first file. Otherwise use this entire set    
    if samples:
        for sample in samples: assert sample in headData["sampleNames"], "Sample {} not in VCF header\n".format(sample)
    else: samples = headData["sampleNames"]
    
    
    ploidyDict = defaultdict(lambda: args.ploidy)
    if args.ploidyFile:
        with open(args.ploidyFile, "rt") as pf: ploidyDict.update(dict([[s[0],int(s[1])] for s in [l.split() for l in pf]]))
    
    ##########################################################################################################################
    
    first_columns = ["#CHROM", "POS"]
    if args.addRefTrack: first_columns.append("REF")
    
    outFile.write(args.outSep.join(first_columns + samples) + "\n")
    
    for vcfSite in parseVcfSites(inFile, headData["mainHeaders"], excludeDuplicates=args.excludeDuplicates, simplifyALT=simplifyALT):
        if (exclude and vcfSite.CHROM in exclude) or (include and vcfSite.CHROM not in include): continue
        if args.minQual and canFloat(vcfSite.QUAL) and float(vcfSite.QUAL) < args.minQual: continue
        if args.maxREFlen and len(vcfSite.REF) > args.maxREFlen: continue
        if args.field is not None: output = vcfSite.getGenoField(args.field,samples=samples, missing=args.missing)
        else:
            output = vcfSite.getGenotypes(gtFilters,asList=True,withPhase=True,samples=samples,missing=args.missing,
                                          mustMatchREFlen=args.skipIndels,keepPartial=args.keepPartial,ploidyDict=ploidyDict,
                                          ploidyMismatchToMissing=args.ploidyMismatchToMissing,expandMulti=args.expandMulti)
        #if we expanded multi-site genotypes, we need to write multiple lines
        if args.expandMulti:
            for x in range(vcfSite.REFlen):
                first_columns = [vcfSite.CHROM, str(vcfSite.POS + x)]
                if args.addRefTrack: first_columns.append(vcfSite.REF[x])
                outFile.write(args.outSep.join(first_columns + [o[x] for o in output]) + "\n")
            continue
        
        first_columns = [vcfSite.CHROM, str(vcfSite.POS)]
        if args.addRefTrack: first_columns.append(vcfSite.REF)
        outFile.write(args.outSep.join(first_columns + output) + "\n")
    
    outFile.close()
