#!/usr/bin/env python

import argparse
import sys
import gzip
import numpy as np
from scipy.stats import ranksums

import genomics

from time import sleep

from threading import Thread

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue



####################################################################################################################################

def which_lowest_test(list_of_arrays, p_threshold = 0.05, noresult=-1):
    n = len(list_of_arrays)
    i = np.argmin([np.nanmean(a) for a in list_of_arrays])
    for j in range(n):
        if i != j:
            result = ranksums(list_of_arrays[i], list_of_arrays[j], alternative = "less")
            if result.pvalue > p_threshold: return noresult
    
    return i

def which_lowest_delta(list_of_arrays, delta_threshold=0, noresult=-1):
    n = len(list_of_arrays)
    means = [np.nanmean(a) for a in list_of_arrays]
    i = np.argmin(means)
    sorted_means = sorted(means)
    
    if sorted_means[1] - sorted_means[0] < delta_threshold: return noresult
    
    return i


'''A function that reads from the window queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def window_wrapper(windowQueue, resultQueue, windType, sampleData, minSites,
                   minData, popNames, ref_pop_indices, p_threshold, delta_threshold=None, noresult=-1, addWindowID=False):
    while True:
        windowNumber,window = windowQueue.get() # retrieve window
        
        if windowNumber == -1:
            resultQueue.put((-1,None,False)) # this is the way of telling everything we're done
            break
        
        if windType == "coordinate" or windType == "predefined":
            scaf,start,end,mid,sites = (window.scaffold, window.limits[0], window.limits[1], window.midPos(),window.seqLen())
        else: scaf,start,end,mid,sites = (window.scaffold, window.firstPos(), window.lastPos(),window.midPos(),window.seqLen())
        
        if sites >= minSites:
            isGood = True
            #make alignment object
            aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = "haplo")
            
            #best match for each ind
            bestMatch = []
            for i in range(len(sampleData.indNames)):
                allPopDists = []
                for popName in popNames:
                    popDists = []
                    for j in ref_pop_indices[popName]:
                        dist = aln.pairDist(i,j)
                        pairSites = np.sum(aln.nanMask[i,:] & aln.nanMask[j,:])
                        popDists.append(dist if pairSites >= minSites else np.nan)
                        
                    allPopDists.append(popDists)
                
                if delta_threshold is not None:
                    bestMatch.append(which_lowest_delta(allPopDists, delta_threshold, noresult))
                else:
                    bestMatch.append(which_lowest_test(allPopDists, p_threshold, noresult))
        
        else:
            isGood = False
            bestMatch = [np.NaN]*len(sampleData.indNames)
        
        results = [] if not addWindowID else [window.ID]
        results += [scaf,start,end,mid,sites] + bestMatch
        resultString = "\t".join([str(x) for x in results])
        resultQueue.put((windowNumber, resultString, isGood))


def sorter(resultQueue, writeQueue, verbose, nWorkerThreads):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
    while True:
        windowNumber,result,isGood = resultQueue.get()
        #check if we're done
        if windowNumber == -1: threadsComplete += 1
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1,None,False))
            break #this is the way of telling everything we're done
        resultsReceived += 1
        if verbose:
            sys.stderr.write("Sorter received window {}\n".format(windowNumber))
        if windowNumber == expect:
            writeQueue.put((windowNumber,result,isGood))
            if verbose:
                sys.stderr.write("Slice {} sent to writer\n".format(windowNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    result,isGood = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,result,isGood))
                    if verbose:
                        sys.stderr.write("Slice {} sent to writer\n".format(expect))
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(windowNumber)] = (result,isGood)



'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out, writeFailedWindows=False):
    global resultsWritten
    global resultsHandled
    while True:
        windowNumber,result,isGood = writeQueue.get()
        #check if we're done
        if windowNumber == -1: break
        if verbose:
            sys.stderr.write("Writer received result {}\n".format(windowNumber))
        if isGood or writeFailedWindows:
            out.write(result + "\n")
            resultsWritten += 1
        resultsHandled += 1


'''loop that checks stats'''
def checkStats():
  while True:
    sleep(10)
    sys.stderr.write("\n{} windows queued, {} results received, {} results written.\n".format(windowsQueued, resultsReceived, resultsWritten))

####################################################################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--windType", help="Type of windows to make", action = "store", choices = ("sites","coordinate","predefined"), default = "coordinate")
    parser.add_argument("-w", "--windSize", help="Window size in bases", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-s", "--stepSize", help="Step size for sliding window", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-m", "--minSites", help="Minumum good sites per window", type=int, action = "store", required = False, metavar="sites", default = 1)
    parser.add_argument("-O", "--overlap", help="Overlap for sites sliding window", type=int, action = "store", required = False, metavar="sites")
    parser.add_argument("-D", "--maxDist", help="Maximum span distance for sites window", type=int, action = "store", required = False)
    parser.add_argument("--windCoords", help="Window coordinates file (scaffold start end)", required = False)
    
    parser.add_argument("--minData", help="Minumum proportion of individuals (or pairs) with >=minSites data", type=float, action = "store", required = False, metavar="prop", default = 0.01)
    
    parser.add_argument("--p_threshold", help="Maximum p-value to assign to a reference populaton", type=float, default=0.05, required = False)
    parser.add_argument("--delta_threshold", help="Alternative method: minimum delta from next best to assign to a reference populaton", type=float, default=None, required = False)
    
    parser.add_argument("-p", "--population", help="Pop name and optionally sample names (separated by commas)",
                        required = False, action='append', nargs="+", metavar=("popName","[samples]"))
    parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)
    parser.add_argument("--samples", help="Samples to include for individual analysis", action = "store", metavar = "sample names")
    
    parser.add_argument("--noresult", help="Value to use when no population is assigned", type=int, default=-1, required = False)
    
    parser.add_argument("-g", "--genoFile", help="Input genotypes file", required = True)
    parser.add_argument("-o", "--outFile", help="Results file", required = False)
    parser.add_argument("--exclude", help="File of scaffolds to exclude", required = False)
    parser.add_argument("--include", help="File of scaffolds to analyse", required = False)
    parser.add_argument("--header", help="Header text if no header in input", action = "store")
    
    parser.add_argument("-T", "--threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("--addWindowID", help="Add window name or number as first column", action="store_true")
    parser.add_argument("--writeFailedWindows", help="Write output even for windows with too few sites.", action="store_true")
    
    
    args = parser.parse_args()
    
    #window parameters
    windType = args.windType
    
    if args.windType == "coordinate":
        assert args.windSize, "Window size must be provided."
        windSize = args.windSize
        stepSize = args.stepSize
        if not stepSize: stepSize = windSize
        assert not args.overlap, "Overlap does not apply to coordinate windows. Use --stepSize instead."
        assert not args.maxDist, "Maximum distance only applies to sites windows."
    
    elif args.windType == "sites":
        assert args.windSize, "Window size (number of sites) must be provided."
        windSize = args.windSize
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
    
    #file info
    
    outFileName = args.outFile
    
    exclude = args.exclude
    include = args.include
    
    #other
    verbose = args.verbose
    
    
    ############## parse samples and populations
    with gzip.open(args.genoFile, "rt") if args.genoFile.endswith(".gz") else open(args.genoFile, "rt") as gf:
        allInds = gf.readline().split()[2:]
    
    popNames = []
    ref_pop_indices = {}
    
    for p in args.population:
        popNames.append(p[0])
        ref_pop_indices[p[0]] = []
        if len(p) > 1:
            for ind in p[1].split(","):
                ref_pop_indices[p[0]].append(allInds.index(ind))
    
    if args.popsFile:
        with open(args.popsFile, "rt") as pf: popDict = dict([ln.split() for ln in pf])
        for ind in popDict.keys():
            try: ref_pop_indices[popDict[ind]].append(allInds.index(ind))
            except: pass
    
    for popName in popNames: assert len(ref_pop_indices[popName]) >= 1, f"Reference population {popName} appears to have no individuals."
    
    ploidyDict = dict(zip(allInds,[1]*len(allInds)))
    
    sampleData = genomics.SampleData(indNames = allInds, ploidyDict = ploidyDict)
    
    ############################################################################################################################################
    
    #open files
    
    genoFile = gzip.open(args.genoFile, "rt") if args.genoFile.endswith(".gz") else open(args.genoFile, "rt")
    
    if args.outFile: outFile = gzip.open(args.outFile, "wt") if args.outFile.endswith(".gz") else open(args.outFile, "wt")
    else: outFile = sys.stdout
    
    if not args.addWindowID: outFile.write("\t".join(["scaffold","start","end","mid","sites"]) + "\t")
    else: outFile.write("\t".join(["windowID","scaffold","start","end","mid","sites"]) + "\t")
    
    ############################################################################################################################################
    
    #stats to output
    
    outFile.write("\t".join(allInds) + "\n")
    
    ##############################################################
    
    #scafs to exclude
    
    if exclude:
        scafsFile = open(exclude, "rU")
        scafsToExclude = [line.rstrip() for line in scafsFile.readlines()]
        sys.stderr.write("{} scaffolds will be excluded.".format(len(scafsToExclude)))
        scafsFile.close()
    else:
        scafsToExclude = None
    
    if include:
        scafsFile = open(include, "rU")
        scafsToInclude = [line.rstrip() for line in scafsFile.readlines()]
        sys.stderr.write("{} scaffolds will be analysed.".format(len(scafsToInclude)))
        scafsFile.close()
    else:
        scafsToInclude = None
    
    
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
    workerThreads = []
    sys.stderr.write("\nStarting {} worker threads\n".format(args.threads))
    for x in range(args.threads):
        workerThread = Process(target=window_wrapper, args = (windowQueue, resultQueue, windType, sampleData, minSites,
                                                             args.minData, popNames, ref_pop_indices,
                                                             args.p_threshold,  args.delta_threshold, args.noresult, args.addWindowID))
        workerThread.daemon = True
        workerThread.start()
        workerThreads.append(workerThread)
    
    
    '''thread for sorting results'''
    sorterThread = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,args.threads,))
    sorterThread.daemon = True
    sorterThread.start()
    
    '''start thread for writing the results'''
    writerThread = Thread(target=writer, args=(writeQueue, outFile, args.writeFailedWindows,))
    writerThread.daemon = True
    writerThread.start()
    
    
    '''start background Thread that will run a loop to check run statistics and print
    We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
    checkerThread = Thread(target=checkStats)
    checkerThread.daemon = True
    checkerThread.start()
    
    
    
    
    ##########################################################
    
    #get windows and analyse
    if windType == "coordinate": windowGenerator = genomics.slidingCoordWindows(genoFile, windSize, stepSize,
                                                                                headerLine = args.header,
                                                                                names = sampleData.indNames,
                                                                                include = scafsToInclude,
                                                                                exclude = scafsToExclude)
    elif windType == "sites": windowGenerator = genomics.slidingSitesWindows(genoFile, windSize, overlap,
                                                                            maxDist, minSites,
                                                                            headerLine = args.header,
                                                                            names = sampleData.indNames,
                                                                            include = scafsToInclude,
                                                                            exclude = scafsToExclude)
    else: windowGenerator = genomics.predefinedCoordWindows(genoFile, windCoords,
                                                            headerLine = args.header,
                                                            names = sampleData.indNames)
    
    
    for window in windowGenerator:
        windowQueue.put((windowsQueued,window))
        windowsQueued += 1
    
    ############################################################################################################################################
    
    #Now we send completion signals to all worker threads
    for x in range(args.threads):
        windowQueue.put((-1,None,)) # -1 tells the threads to break
    
    sys.stderr.write("\nWaiting for all threads to finish\n".format(args.threads))
    for x in range(len(workerThreads)):
        workerThreads[x].join()
    
    sorterThread.join()
    writerThread.join()
    
    genoFile.close()
    outFile.close()
    
    sys.stderr.write(str(windowsQueued) + " windows were tested.\n")
    sys.stderr.write(str(resultsWritten) + " results were written.\n")
    
    sys.stderr.write("\nDone.\n")
    
    sys.exit()



