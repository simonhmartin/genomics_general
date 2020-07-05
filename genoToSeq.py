#script converts genotype format files to sequence alignments (fasta or phylip)
import gzip, argparse, sys
import genomics


#############################

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genoFile", help="Input geno file", action = "store", required = False)
parser.add_argument("-s", "--seqFile", help="Output sequence file", action = "store", required = False)
parser.add_argument("-f", "--format", help="Sequence file format", action = "store", required = False, choices = ("phylip","fasta",), default = "fasta")
parser.add_argument("-M", "--mode", help="Output mode for different contigs", action = "store",
                    choices = ("cat", "windows","contigs"), default = "cat")
parser.add_argument("-S", "--samples", help="Name of sample(s)", action = "store", required = False)
parser.add_argument("--NtoGap", help="Convert 'N' or 'n' to '-'", action = "store_true")
parser.add_argument("--seqNameFormat", help="Format for sequence names", action = "store", required = False,
                    choices = ("sample", "contig", "sample_contig", "contig_position", "sample_contig_position"),  default = "sample")
parser.add_argument("--splitPhased", help="Split phased genotypes into two (or more, see --ploidy) sequences per sample", action = "store_true")
parser.add_argument("--ploidy", help="Ploidy for each individual. Only necessary when splitting phased sequences", action = "store", nargs="+", type=int, default = [2])
parser.add_argument("--separateFiles", help="Output windows or contigs as separate files", action = "store_true")
parser.add_argument("--gzip", help="gzip output file(s)", action = "store_true")

#arguments for sliding windows
parser.add_argument("--windType", help="FOR WINDOWS: type of windows to make", action = "store", choices = ("sites","coordinate"), default = "sites")
parser.add_argument("--windSize", help="FOR WINDOWS: Window size in bases", type=int, action = "store")
parser.add_argument("--minSites", help="FOR WINDOWS: Minumum sites per window", type=int, action = "store")
parser.add_argument("--stepSize", help="FOR WINDOWS: Step size for coordinate sliding window", type=int, action = "store")
parser.add_argument("--overlap", help="FOR WINDOWS: Overlap for sites sliding window", type=int, action = "store")
parser.add_argument("--maxDist", help="FOR WINDOWS: Maximum span distance for sites window", type=int, action = "store")




args = parser.parse_args()

#open input
if args.genoFile:
    if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "rt")
    else: genoFile = open(args.genoFile, "rt")
else: genoFile = sys.stdin

#open output (unless doing separate uoutputs for all contigs
if not args.separateFiles:
    if args.seqFile:
        if args.seqFile[-3:] == ".gz": seqFile = gzip.open(args.seqFile, "wt")
        elif args.gzip: seqFile = gzip.open(args.seqFile+".gz", "wt")
        else: seqFile = open(args.seqFile, "wt")
    else: seqFile = sys.stdout

#############################

samples = args.samples.split(",") if args.samples else None

# if cating all contigs, just parse file and write
if args.mode == "cat":
    #read file into window like object
    window = genomics.parseGenoFile(genoFile, names=samples, splitPhased=args.splitPhased, ploidy=args.ploidy)
    #write
    seqDict = window.seqDict()
    seqFile.write(genomics.makeAlnString(window.names,[seqDict[name] for name in window.names],outFormat = args.format, NtoGap=args.NtoGap))
    genoFile.close()
    seqFile.close()
    exit()

if args.mode == "windows" or args.mode == "contigs":
    if args.mode == "windows":
        windType = args.windType
        windSize = args.windSize
        minSites = args.minSites
        stepSize = args.stepSize
        overlap = args.overlap
        maxDist = args.maxDist
    else:
        #to get contigs, we just use very lare non-overlapping coordinate windows
        windType = "coordinate"
        windSize = 1e7
        minSites = 1
        stepSize = 1e7
    #initialise window generator
    if windType == "coordinate": windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize,args.samples,
                                                                                splitPhased=args.splitPhased, ploidy=args.ploidy)
    else: windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap, maxDist, minSites, args.samples,
                                                         splitPhased=args.splitPhased, ploidy=args.ploidy)
    
    for window in windowGenerator:
        
        posString = str(window.firstPos()) + "_" + str(window.lastPos())
        
        #open file if necessary
        if args.separateFiles:
            if args.format == "fasta": ext = ".fa"
            else: ext = ".phy"
            seqFileName = args.seqFile + "." + window.scaffold
            if args.mode == "windows": seqFileName += "_" + posString
            seqFileName += ext
            if args.gzip:
                seqFileName += ".gz"
                seqFile = gzip.open(seqFileName, "wt")
            else: seqFile = open(seqFileName, "wt")
        
        #change names if necessary
        
        if args.seqNameFormat == "contig":
            seqNames = [window.scaffold]*window.n
        elif args.seqNameFormat == "sample_contig":
            seqNames = [n+"_"+window.scaffold for n in window.names]
        elif args.seqNameFormat == "contig_position":
            seqNames = [window.scaffold + "_" + posString]*window.n
        elif args.seqNameFormat == "sample_contig_position":
            seqNames = [n+"_"+window.scaffold+"_"+posString for n in window.names]
        else: seqNames = window.names
        
        #write
        seqDict = window.seqDict()
        seqFile.write(genomics.makeAlnString(seqNames,[seqDict[name] for name in seqNames],outFormat = args.format,NtoGap=args.NtoGap))
        
        if args.separateFiles: seqFile.close()

seqFile.close()
genoFile.close()

exit()
