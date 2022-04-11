#!/usr/bin/env python

import argparse, itertools, sys, gzip, numpy

import genomics

from time import sleep

### parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Input geno file", action = "store")
parser.add_argument("-f", "--genoFormat", help="Data format of input", action = "store", choices = ("phased","diplo","alleles"), default = "phased")
parser.add_argument("-o", "--outfile", help="Output csv file", action = "store")
parser.add_argument("-s", "--samples", help="Specify samples to include", action = "store", metavar = "comma separated names")
parser.add_argument("--maxAlleles", help="Maximum number of alleles per site", type=int, action = "store", default = 2, choices = [2,3,4])
parser.add_argument("--includeNull", help="Include null genotypes", action = "store_true")
parser.add_argument("--maxSites", help="Maximum number of sites to count", type=int, action = "store")

parser.add_argument("--test", help="Test run", action = "store_true")



args = parser.parse_args()

samples = args.samples.split(",") if args.samples else None

test = args.test

### open files

if args.infile:
    In = gzip.open(args.infile, "rt") if args.infile.endswith(".gz") else open(args.infile, "rt")    
else:
    In = sys.stdin

if args.outfile:
    Out = gzip.open(args.outfile, "wt") if args.outfile.endswith(".gz") else open(args.outfile, "wt")
else:
    Out = sys.stdout


### parse vcf header

reader = genomics.GenoFileReader(In)

if samples:
    for sample in samples:
        assert sample in reader.names, "Specified sample name not in VCF header."
else:
    samples = reader.names

nSamples = len(samples)


sys.stderr.write("\n {} samples will be considered.".format(nSamples))


### generate code combinations

elements = [str(x) for x in range(args.maxAlleles)]

if args.includeNull:
    elements += ["N"]

genotypes = ["".join(x) for x in list(itertools.combinations_with_replacement(elements, 2))]


sys.stderr.write("\nThe following genotypes will be considered:\n")
sys.stderr.write(" ".join(genotypes))

### check that we aren't going to generate too many combinations

nPatterns = len(genotypes)**nSamples

sys.stderr.write("\nThis corresponds to {} unique patterns.\n".format(nPatterns))
assert nPatterns <= 1000000, "Trying to evaluate this many patterns will use too much memory."


patterns = list(itertools.product(genotypes, repeat = nSamples))

#make a dictionary of patterns, this will be MUCH faster to search

patDict = dict(zip(patterns, range(nPatterns)))


#and variable to record counts of each

counts = numpy.zeros(nPatterns)


### count patterns

done = 0

for siteData in reader.siteBySite():
    if test:
        sys.stderr.write("{},{},{}\n".format(siteData["scaffold"], siteData["positin"], siteData["GTs"]))
        sleep(0.05)
    site = genomics.GenomeSite(genotypes=[siteData["GTs"][name] for name in samples], sampleNames=samples,
                               genoFormat=args.genoFormat)
    pattern = tuple([gt[::2] for gt in site.asList(mode = "coded")])
    if test: print(",".join(pattern), file=sys.stderr)
    if pattern in patDict: counts[patDict[pattern]] += 1
    
    done += 1
    if args.maxSites and done >= args.maxSites:
        break

In.close()


### wriute output

Out.write(",".join(samples) + ",count\n")

for n in range(nPatterns):
    for geno in patterns[n]:
        Out.write(geno + ",")
    Out.write(str(counts[n]) + "\n")

Out.close()

    
    
