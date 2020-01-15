#!/usr/bin/env python

import argparse, sys, os, gzip, random, tempfile
import numpy as np

from threading import Thread

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue

from time import sleep

#add parent directory to path in case genomics.py is there instead of on the path
sys.path.insert(1, os.path.realpath(os.path.pardir))
import genomics


####################################################################################################################################

#command that actually runs phyml and gets the bestTree 
def phymlTree(seqArray, seqNames, model, opt, phyml, prefix = "", tmpDir = None, test = False, log="/dev/null"):
    #write file
    tempAln = tempfile.NamedTemporaryFile(mode="w",prefix=prefix,suffix=".phy",dir=tmpDir,delete=False)
    localName = tempAln.name.rsplit("/",1)[1]
    with tempAln as tA: tA.write(genomics.makeAlnString(seqNames,seqArray))
    phymlCommand = " ".join([phyml,"--input", tempAln.name,"--model", model, "-o", opt, "-b 0", ">>", log])
    if test: sys.stderr.write("phyml command:\n" + phymlCommand + "\n")
    os.system(phymlCommand)
    #try retrieve the result  
    try:
        with open(tempAln.name + "_phyml_tree.txt", "rt") as treeFile: tree = treeFile.readline().strip()
    except:
        try:
            with open(tempAln.name + "_phyml_tree", "rt") as treeFile: tree = treeFile.readline().strip()
        except:
            if verbose:
                sys.stderr.write("Tree not found at " + tempAln.name + "_phyml_tree.txt\n")
            tree = "NA"
    try:
        with open(tempAln.name + "_phyml_stats.txt", "rt") as statsFile:
            stats = statsFile.read().split()
            lnL = stats[stats.index("Log-likelihood:")+1]
    except:
        try:
            with open(tempAln.name + "_phyml_stats", "rt") as statsFile:
                stats = statsFile.read().split()
                lnL = stats[stats.index("Log-likelihood:")+1]
        except:
            lnL = "NA"
    #remove files
    if not test:
        for f in [f for f in os.listdir(tmpDir) if localName in f]:
            os.remove(tmpDir + "/" + f)
    return (tree,lnL,)


#command that actually runs phyml and gets the bestTree 
def phymlCrossVal(seqArray0, seqArray1, indNames, model, opt, phyml, prefix = "",tmpDir=None, test = False, log="/dev/null"):
    #write file
    tempAln0 = tempfile.NamedTemporaryFile(mode="w",prefix=prefix,suffix=".0.phy",dir=tmpDir,delete=False)
    tempAln1 = tempfile.NamedTemporaryFile(mode="w",prefix=prefix,suffix=".1.phy",dir=tmpDir,delete=False)
    localName0 = tempAln0.name.rsplit("/",1)[1]
    localName1 = tempAln1.name.rsplit("/",1)[1]
    with tempAln0 as tempAln0: tempAln0.write(genomics.makeAlnString(seqNames,seqArray0))
    with tempAln1 as tempAln1: tempAln1.write(genomics.makeAlnString(seqNames,seqArray1))
    #first way validation
    #tree
    phymlCommand = " ".join([phyml,"--input", tempAln0.name,"--model", model, "-o", opt, ">>", log])
    os.system(phymlCommand)
    #validation
    phymlCommand = " ".join([phyml,"--input", tempAln1.name,"--model", model, "-o", "n", "-u", tempAln0.name + "_phyml_tree.txt", ">>", log])
    os.system(phymlCommand)
    #retrieve
    try:
        with open(tempAln1.name + "_phyml_stats.txt", "rt") as statsFile:
            stats = statsFile.read().split()
            lnL1 = float(stats[stats.index("Log-likelihood:")+1])
    except: lnL1 = np.NaN
    #second way validation
    #tree
    phymlCommand = " ".join([phyml,"--input", tempAln1.name,"--model", model, "-o", opt, ">>", log])
    os.system(phymlCommand)
    #validation
    phymlCommand = " ".join([phyml,"--input", tempAln0.name,"--model", model, "-o", "n", "-u", tempAln1.name + "_phyml_tree.txt", ">>", log])
    os.system(phymlCommand)
    #retrieve
    try:
        with open(tempAln0.name + "_phyml_stats.txt", "rt") as statsFile:
            stats = statsFile.read().split()
            lnL0 = float(stats[stats.index("Log-likelihood:")+1])
    except: lnL0 = np.NaN
    #remove files
    if not test:
        for f in [f for f in os.listdir(tmpDir) if localName0 in f or localName1 in f]:
            os.remove(tmpDir + "/" + f)
    return str(lnL0+lnL1)


'''A function that reads from the window queue, calls sume other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using'''
def phyml_wrapper(windowQueue, resultQueue, windType, model, opt, outgroup, phyml, minSites, minPerInd, minSNPs=None,
                  maxLDphase=False, bootstraps=0, crossVal=False, test = False):
    while True:
        windowNumber,window = windowQueue.get()
        Nsites = window.seqLen()
        if test or verbose: sys.stderr.write("Window {} received for analysis, length: {}\n".format(windowNumber,Nsites))
        if windType == "coordinate" or windType == "predefined": scaf,start,end,mid = (window.scaffold, window.limits[0], window.limits[1], window.midPos())
        else: scaf,start,end,mid = (window.scaffold, window.firstPos(), window.lastPos(), window.midPos())
        prefix = scaf + "_" + str(start) + "_" + str(end) + "_"
        if Nsites >= minSites:
            
            aln = genomics.genoToAlignment(window.seqDict(), genoFormat = "phased")
            if len(outgroup) >= 1:
                for seqName in aln.names:
                    if seqName in outgroup: seqName +="*"
            
            sitesPerInd = aln.seqNonNan()
            if min(sitesPerInd) >= minPerInd and (minSNPs is None or
                                                  len(aln.varSites(indices=np.array([i for i in range(aln.N) if aln.sampleNames[i] not in outgroup]))) >= minSNPs):
                if maxLDphase: aln = genomics.maxLDphase(aln)
                #if enough sites get tree
                tree,lnL = phymlTree(aln.array,aln.names,model,opt,phyml,prefix,tmpDir=tmpDir, test = test, log = log)
                bsTrees = []
                for b in range(bootstraps):
                    #get bootstrap trees if necessary
                    positions = np.random.choice(range(Nsites), Nsites, replace=True)
                    newArr = aln.array[:,positions]
                    bsTree,bslnL = phymlTree(newArr,aln.names,model,opt,phyml,prefix + str(b) + "_",tmpDir=tmpDir, test = test, log = log)
                    bsTrees.append(bsTree)
                trees = [tree] + bsTrees
                if crossVal:
                    arr0 = aln.arr[:,range(int(round(Nsites/2)))]
                    arr1 = aln.arr[:,range(int(round(Nsites/2)), Nsites)]
                    cvlnL = phymlCrossVal(arr0,arr1,indNames,model,opt,phyml,prefix,tmpDir=tmpDir, test = test, log = log)
            else:
                trees = ["NA"] + ["NA"]*bootstraps
                lnL = cvlnL = "NA"
        else:
            trees = ["NA"] + ["NA"]*bootstraps
            lnL = cvlnL = "NA"
                
        data = [window.scaffold, str(start), str(end), str(mid), str(Nsites), str(lnL)]
        if crossVal: data.append(cvlnL)
        
        output = ["\t".join(data)] + trees
        
        resultQueue.put((windowNumber, tuple(output),))


'''a function that watches the result queue and sorts results.'''
def sorter(resultQueue, writeQueue, verbose):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    while True:
        resNumber,result = resultQueue.get()
        resultsReceived += 1
        if verbose: sys.stderr.write("Sorter received result " + str(resNumber))
        if resNumber == expect:
            writeQueue.put((resNumber,result,))
            if verbose: sys.stderr.write("Result {} sent to writer".format(resNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try: result = sortBuffer.pop(str(expect))
                except: break
                #if we get here we've found the one we want in the buffer
                writeQueue.put((expect,result))
                if verbose: sys.stderr.write("Result {} sent to writer".format(expect))
                expect +=1
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(resNumber)] = result


'''a writer function that writes the sorted result.'''
def writer(writeQueue, outs):
    global resultsWritten
    global resultsHandled
    while True:
        resNumber,result = writeQueue.get()
        if verbose: sys.stderr.write("Writer received result {}\n".format(resNumber))
        for x in range(len(outs)): outs[x].write(result[x] + "\n")
        resultsWritten += 1
        resultsHandled += 1

'''loop that checks stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("\n{} windows queued | {} results received | {} results written.\n".format(windowsQueued, resultsReceived, resultsWritten))

####################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("--windType", help="Type of windows to make", action = "store", choices = ("sites","coordinate","predefined"), default = "coordinate")

parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-M", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-Mi", "--minPerInd", help="Minumum good sites per individual", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-Ms", "--minSNPs", help="Minumum number of SNPs", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-S", "--stepSize", help="Step size for coordinate sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)
parser.add_argument("--windCoords", help="Window coordinates file (scaffold start end)", required = False)

parser.add_argument("-g", "--genoFile", help="Input genotypes file")
parser.add_argument("-p", "--prefix", help="Prefix for output files", required = True)

parser.add_argument("--exclude", help="List of scaffolds to exclude (comma separated)", required = False)
parser.add_argument("--include", help="List of scaffolds to analyse (comma separated)", required = False)
parser.add_argument("--excludeFile", help="File of scaffolds to exclude", required = False)
parser.add_argument("--includeFile", help="File of scaffolds to analyse", required = False)

parser.add_argument("--individuals", help="Individuals to include, separated by comma", action = "store", metavar="ind1,ind2,ind3...")

parser.add_argument("--indFile", help="File of individuals to include, one per line", action = "store")

parser.add_argument("--maxLDphase", help="Do crude phase improvement based on LD", action = "store_true")

parser.add_argument("--outgroup", help="Outgroup individuals, separated by comma", action = "store", metavar="ind1,ind2,ind3...")

parser.add_argument("--phyml", help="path to phyml executable", action = "store", metavar="path/to/phyml", default="phyml")
parser.add_argument("--model", help="phyml model", action = "store", default="GTR")
parser.add_argument("--optimise", help="parameters to optimise - see phyml manual", action = "store", choices = ("tlr","tl","tr","lr","t","l","r","n"), default="n")

parser.add_argument("--bootstraps", help="number of bootstrap resamplings to do", type=int, action = "store", default=0)
parser.add_argument("--crossVal", help="do cross validation and report likelihood", action = "store_true")

parser.add_argument("--tmp", help="Location for temporary phyml files", action = "store", metavar="path/to/tmp")
parser.add_argument("--log", help="phyml log file, if you want one.", action = "store", metavar="path/to/log", default="/dev/null")

parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
parser.add_argument("--test", help="run test", action="store_true")
parser.add_argument("--verbose", help="Verbose output", action="store_true")


args = parser.parse_args()

verbose = args.verbose
test = args.test

windType = args.windType

windSize = args.windSize

if args.windType == "coordinate":
    assert args.windSize, "Window size must be provided."
    stepSize = args.stepSize
    if not stepSize: stepSize = windSize
    assert not args.overlap, "Overlap noes not apply to coordinate windows. Use --stepSize instead."
    assert not args.maxDist, "Maximum distance only applies to sites windows."

elif args.windType == "sites":
    assert args.windSize, "Window size (number of sites) must be provided."
    overlap = args.overlap
    if not overlap: overlap = 0
    maxDist = args.maxDist
    if not maxDist: maxDist = np.inf
    assert not args.stepSize, "Step size only applies to coordinate windows. Use --overlap instead."
else:
    assert args.windCoords, "Please provide a file of window coordinates."
    assert not args.overlap, "Overlap does not apply for predefined windows."
    assert not args.maxDist, "Maximum does not apply for predefined windows."
    assert not args.stepSize,"Step size does not apply for predefined windows."
    assert not args.include,"You cannot only include specific scaffolds if using predefined windows."
    assert not args.exclude,"You cannot exclude specific scaffolds if using predefined windows."
    with open(args.windCoords,"rt") as wc: windCoords = tuple([(x,int(y),int(z),) for x,y,z in [line.split()[:3] for line in wc]])

minSites = args.minSites
if not minSites: minSites = windSize

minPerInd = args.minPerInd if args.minPerInd else minSites


if args.individuals: indNames = args.individuals.split(",")
elif args.indFile:
    with open(args.indFile, "rt") as indFile:
        indNames = [name.strip() for name in indFile.readlines()]
else: indNames = None

if args.outgroup:
    outgroup = args.outgroup.split(",")
    if test or verbose: sys.stderr.write("outgroups:" + " ".join(outgroup) + "\n")
else: outgroup = []

log = args.log

threads = args.threads

phyml = args.phyml
model = args.model
opt = args.optimise

bootstraps = args.bootstraps

############################################################################################################################################

#open files

if args.genoFile: genoFile = gzip.open(args.genoFile, "rt") if args.genoFile.endswith(".gz") else open(args.genoFile, "rt")
else: genoFile = sys.stdin

dataFile = open(args.prefix + ".data.tsv", "wt")

outHeads = ["scaffold","start","end","mid","sites","lnL"]
if args.crossVal: outHeads.append("cv_lnL")

dataFile.write("\t".join(outHeads) + "\n")

treesFile = gzip.open(args.prefix + ".trees.gz", "wt")

outs = [dataFile, treesFile]

for b in range(bootstraps): outs.append(gzip.open(args.prefix + ".BS" + str(b) + ".trees.gz", "wt"))

#tmp dir for phyml work

tmpDir = tempfile.mkdtemp(prefix="phyml_tmp", dir=args.tmp)
if test or verbose: sys.stderr.write("\nTemporary Phyml files will be stored in {}\n".format(tmpDir))

############################################################################################################################################


#scafs to exclude

if args.exclude:
    scafsToExclude = args.exclude.split(",")
    sys.stderr.write("{} scaffolds will be excluded.".format(len(scafsToExclude)))
elif args.excludeFile:
    with open(args.excludeFile, "rt") as scafsFile: scafsToExclude = [line.rstrip() for line in scafsFile]
    sys.stderr.write("{} scaffolds will be excluded.".format(len(scafsToExclude)))
else: scafsToExclude = None

if args.include:
    scafsToInclude = args.include.split(",")
    sys.stderr.write("{} scaffolds will be analysed.".format(len(scafsToInclude)))
elif args.includeFile:
    with open(args.includeFile, "rt") as scafsFile: scafsToInclude = [line.rstrip() for line in scafsFile]
    sys.stderr.write("{} scaffolds will be analysed.".format(len(scafsToInclude)))
else: scafsToInclude = None


##########################################################################################################

#counting stat that will let keep track of how far we are
windowsQueued = 0
resultsReceived = 0
resultsWritten = 0
resultsHandled = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
windowQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''

for x in range(threads):
    worker = Process(target=phyml_wrapper, args = (windowQueue, resultQueue, windType,
                                                   model, opt, outgroup, phyml, minSites, minPerInd,
                                                   args.minSNPs, args.maxLDphase, bootstraps, args.crossVal, test,))
    worker.daemon = True
    worker.start()
    if test or verbose: sys.stderr.write("started worker {}\n".format(x))
    

'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, outs,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()




##########################################################

#get windows and analyse
if windType == "coordinate":
    windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize, indNames, include = scafsToInclude, exclude = scafsToExclude)
elif windType == "sites":
    windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap, maxDist, minSites, indNames, include = scafsToInclude, exclude = scafsToExclude)
else:
    windowGenerator = genomics.predefinedCoordWindows(genoFile, windCoords, indNames)

for window in windowGenerator:
    #simpleque has no max, so to make sure we haven't gotten ahead of ourselves, we compare windowsQueued to resultsReceived
    while windowsQueued - resultsReceived >= 50:
        sleep(10)
        if test or verbose: sys.stderr.write("Waiting for queue to clear...\n")
    
    if test or verbose:
        if test: sleep(0.5)
        sys.stderr.write("Sending window {} to queue. Length: {}\n".format(windowsQueued, window.seqLen()))
    
    windowQueue.put((windowsQueued,window))
    windowsQueued += 1
    if test and windowsQueued == 10: break

############################################################################################################################################

sys.stderr.write("\nWriting final results...\n")
while resultsHandled < windowsQueued:
  sleep(1)

sleep(5)

for o in outs: o.close()

sys.stderr.write(str(windowsQueued) + " windows were tested.\n")
sys.stderr.write(str(resultsWritten) + " results were written.\n")

if not test: os.rmdir(tmpDir)

sys.exit()

