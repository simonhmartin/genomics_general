import sys, argparse, gzip

#############

def openForReading(fileName):
    return gzip.open(fileName, "rt") if fileName.endswith(".gz") else open(fileName, "rt")

def openForWriting(fileName):
    return gzip.open(fileName, "w") if fileName.endswith(".gz") else open(fileName, "w")

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", help="Input file", action = "append", required = True)
parser.add_argument("-f", "--fai", help="Reference fasta index file", action = "store", required = True)
parser.add_argument("-o", "--outputFile", help="Output file", action = "store")
parser.add_argument("--method", help="How to merge", action = "store", choices = ("intersect","union","all"), default = "intersect")
parser.add_argument("--unionMin", help="Minimum files represented for menthod union", action = "store", type=int, default=1)
parser.add_argument("--mustIncludeFirst", help="The first n files MUST be present", action = "store", type=int, default=0)
parser.add_argument("--outSep", help="Output file separator", action = "store", default = "\t")
parser.add_argument("--missing", help="Missing genotype for method union or all", action = "store", default = "N")
parser.add_argument("--outputOnly", help="Which output files to include in output", action = "store", type=int, nargs="+")
parser.add_argument('--verbose', help="Verbose output", action='store_true')

args = parser.parse_args()

#create a list of file handles
files = [openForReading(f) for f in args.inputFile]
nFiles = len(files)

out = openForWriting(args.outputFile) if args.outputFile else sys.stdout

#index of inpuit files to include in output (by default this is all)
outputIdx = [i-1 for i in args.outputOnly] if args.outputOnly else range(nFiles) 

#get scaffold list and lengths
with open(args.fai, "rt") as fai: scafLens = [(s,int(l)) for s,l in [ln.split()[:2] for ln in fai]]
scafs = [x[0] for x in scafLens]
scafLens = dict(scafLens)

linesWritten = 0

headers = [file.readline().split() for file in files]

#make dummy lines for each file if necessary
dummyGenos = [[args.missing]*(len(h)-2) for h in headers]

#check that minumum for union is at least as large as the mustIncludeFirst value
unionMin = args.unionMin
if unionMin < args.mustIncludeFirst: unionMin = args.mustIncludeFirst

#write headers
out.write(args.outSep.join([args.outSep.join(headers[0][0:2]),
                            args.outSep.join([args.outSep.join(headers[x][2:]) for x in outputIdx])]) + "\n")

#now start reading all files
objectsList = [file.readline().split() for file in files]

for scaf in scafs:
    sys.stderr.write("Merging {}...\n".format(scaf))
    for site in range(1,scafLens[scaf]+1):
        site = str(site)
        filesRepresented = 0
        outObjects = [scaf, site]
        fail = False
        for x in range(nFiles):
            if len(objectsList[x]) >= 2 and scaf == objectsList[x][0] and site == objectsList[x][1]:
            #if [scaf,site] == objectsList[x][:2]:
                #if its a match, add the objects to the output and read in the next line for that file
                if x in outputIdx: outObjects += objectsList[x][2:]
                objectsList[x] = files[x].readline().split()
                filesRepresented += 1
            else:
                #if not a match, add Ns for this file (if necessary), and dont read next line
                #but can't break because we need to move to next for those that do have the current line
                #but we can record that it's a fail and move to the next line
                if args.method == "intersect" or x < args.mustIncludeFirst:
                    fail = True
                    continue
                if x in outputIdx: outObjects += dummyGenos[x]
        
        #so now we've created the output, but need to decide if we can write it
        if args.verbose: sys.stderr.write("{} {}: {} files represented.\n".format(scaf,site,filesRepresented))
        if not fail and (args.method == "all" or (args.method == "union" and filesRepresented >= unionMin) or (args.method == "intersect" and filesRepresented == nFiles)):
            if args.verbose: sys.stderr.write("Writing line.\n")
            out.write(args.outSep.join(outObjects) + "\n")
            linesWritten += 1
            if linesWritten % 100000 == 0:
                sys.stderr.write("{} lines written to output...\n".format(linesWritten))
        #and thats it. Move on to the next site in the genome

for f in files:
  f.close
