#!/usr/bin/env python

import argparse, sys, gzip, random
import numpy as np

from threading import Thread

from multiprocessing import Process

if sys.version_info>=(3,0):
    from multiprocessing import SimpleQueue
else:
    from multiprocessing.queues import SimpleQueue

import genomics

from time import sleep


#######################################################################################################################

#function to generate slices of a file. N is the sizehint, in bytes I guess, so 1M typically translates fto a few thosand lines
def fileSlicer(f, N=1000000):
    while True:
        fileSlice = f.readlines(N)
        if len(fileSlice) != 0: yield fileSlice
        else: break


'''A function that reads from the input queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def freqs_wrapper(inQueue, resultQueue, headerLine, genoFormat, sampleData, target, minData, asCounts, threshold, keepNanLines = False):
    while True:
        
        sliceNumber,fileSlice = inQueue.get() # retrieve slice
        
        if sliceNumber == -1:
            resultQueue.put((-1,None,)) # this is the way of telling everything we're done
            break

        
        window = genomics.parseGenoFile(fileSlice, headerLine, names=sampleData.indNames)
        
        #make alignment objects
        aln = genomics.genoToAlignment(window.seqDict(), sampleData, genoFormat = genoFormat)
        popAlns = dict([(popName, aln.subset(groups=[popName])) for popName in sampleData.popNames])
        #this above replaced this below, as it should be faster
        #popAlns = dict(zip(sampleData.popNames, [aln.subset(groups=[pop]) for pop in sampleData.popNames]))
        
        #if there is no target, fetch all base counts
        
        if not target:
            popFreqs = []
            for pop in sampleData.popNames:
                goodData = popAlns[pop].siteNonNan() >= minData
                sites = np.where(goodData)[0]
                baseFreqs = popAlns[pop].siteFreqs(asCounts=asCounts)
                popFreqs.append([",".join(row) for row in baseFreqs.astype(str)])
            
            allFreqs = np.column_stack(popFreqs)
            
        else:
            #otherwise define the target base at each site
            if target == "derived":
                #use last pop as outgroup
                outgroup = sampleData.popNames[-1]
                inAln = aln.subset(groups = sampleData.popNames[:-1])
                baseColumns = np.array([genomics.derivedAllele(inAln.numArray[:,i][inAln.nanMask[:,i]],
                                                            popAlns[outgroup].numArray[:,i][popAlns[outgroup].nanMask[:,i]],
                                                            numeric=True)
                                        for i in range(aln.l)]).reshape([aln.l,1])
                
            else:
                #otherwise get minor allele.
                baseColumns = np.array([genomics.minorAllele(aln.numArray[:,i][aln.nanMask[:,i]]) for i in xrange(aln.l)]).reshape([aln.l,1])
            
            goodSites = np.apply_along_axis(lambda x: ~np.any(np.isnan(x)),1,baseColumns)
            
            #get freqs per pop
            popFreqs = []
            for pop in sampleData.popNames:
                #first find sites with sufficient data
                goodData = popAlns[pop].siteNonNan() >= minData
                sites = np.where(goodSites & goodData)[0]
                baseFreqs = popAlns[pop].siteFreqs(sites, asCounts=asCounts)
                popColumns = baseColumns[sites,:].astype(int)
                popRows = np.repeat(np.arange(len(sites))[:,np.newaxis],popColumns.shape[1], axis = 1)
                targetFreqs =  np.zeros([aln.l, popColumns.shape[1]], dtype=int if asCounts else float)
                if not asCounts: targetFreqs.fill(np.nan)
                if len(sites) >= 1: targetFreqs[sites,:] = baseFreqs[popRows,popColumns]
                popFreqs.append(np.around(targetFreqs, 4))
            
            allFreqs = np.hstack(popFreqs)
            
            if threshold and not asCounts:
                allFreqs[allFreqs >= threshold] = 1
                allFreqs[allFreqs < threshold] = 0
        
        #fetch scaffold and position
        scafPos = np.array([line.split(None, 2)[:2] for line in fileSlice], dtype="str")
        
        if not keepNanLines:
            if not asCounts:
                outSites = np.where(~np.apply_along_axis(np.all, 1, np.isnan(allFreqs)))[0]
            else: outSites = np.where(~np.apply_along_axis(np.all, 1, allFreqs==0))[0]
        else: outSites = range(aln.l)
                
        outArray = np.column_stack((scafPos[outSites,:],
                                    allFreqs[outSites,:].astype(str),))
        
        resultStrings = ["\t".join(row) for row in outArray]
        
        resultQueue.put((sliceNumber, resultStrings,))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(doneQueue, writeQueue, verbose, nWorkerThreads):
    global resultsReceived
    sortBuffer = {}
    expect = 0
    threadsComplete = 0 #this will keep track of the worker threads and once they're all done this thread will break
    while True:
        sliceNumber, results = doneQueue.get()
        #check if we're done
        if sliceNumber == -1: threadsComplete += 1
        if threadsComplete == nWorkerThreads:
            writeQueue.put((-1,None,))
            break #this is the way of telling everything we're done
        resultsReceived += 1
        if verbose:
            sys.stderr.write("Sorter received slice {}\n".format(sliceNumber))
        if sliceNumber == expect:
            writeQueue.put((sliceNumber,results))
            if verbose:
                sys.stderr.write("Slice {} sent to writer\n".format(sliceNumber))
            expect +=1
            #now check buffer for further results
            while True:
                try:
                    results = sortBuffer.pop(str(expect))
                    writeQueue.put((expect,results))
                    if verbose:
                        sys.stderr.write("Slice {} sent to writer\n".format(expect))
                    expect +=1
                except:
                    break
        else:
            #otherwise this line is ahead of us, so add to buffer dictionary
            sortBuffer[str(sliceNumber)] = results



'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out, verbose):
    global resultsWritten
    global linesWritten
    while True:
        sliceNumber, results = writeQueue.get()
        #check if we're done
        if sliceNumber == -1: break
        if verbose:
            sys.stderr.write("Writer received slice {}\n".format(sliceNumber))
        for outLine in results:
            out.write(outLine + "\n")
            linesWritten += 1
        resultsWritten += 1

'''loop that checks line stats'''
def checkStats():
    while True:
        sleep(10)
        sys.stderr.write("{} slices queued | {} slices analysed | {} slices written | {} lines written\n".format(slicesQueued,resultsReceived,resultsWritten,linesWritten))

def lineReader(fileObj):
    line = fileObj.readline()
    while len(line) >= 1:
        yield line
        line = fileObj.readline()

#########################################################################################################################

if __name__ == '__main__':

    ### parse arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genoFile", help="Input vcf file", action = "store")
    parser.add_argument("-o", "--outFile", help="Output csv file", action = "store")
    parser.add_argument("-f", "--genoFormat", help="Data format for output", action = "store", choices = ("phased","diplo","alleles"), default = "phased")

    #populations
    parser.add_argument("-p", "--population", help="Pop name and optionally sample names (separated by commas)",
                        required = True, action='append', nargs="+", metavar=("popName","[samples]"))
    parser.add_argument("--popsFile", help="Optional file of sample names and populations", action = "store", required = False)

    #Frequency for all bases or a single base
    parser.add_argument("--target", help="All or single base frequency (derived assumes last pop is outgroup)", choices = ("minor","derived"),
                        action = "store", default=None)

    parser.add_argument("--asCounts", help="Return frequencies as counts", action='store_true')

    #define ploidy if not 2
    parser.add_argument("--ploidy", help="Ploidy for each sample", action = "store", type=int, nargs="+")
    parser.add_argument("--ploidyFile", help="File with samples names and ploidy as columns", action = "store")
    parser.add_argument("--haploid", help="Alternative way to specify samples that are not diploid", action = "store", nargs="+")

    #optional missing data argument
    parser.add_argument("--minData", help="Minimum proportion of non-missing data per population", type=float, action = "store", default = 0, metavar = "proportion")

    #threshold value for rounding to 0 or 1 (only for very specific applicatons)
    parser.add_argument("--threshold", help="Threshold value for rounding to 0 or 1", type=float, action = "store", metavar = "proportion")

    parser.add_argument("--keepNanLines", help="Output lines with no information", action='store_true')


    #other
    parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
    parser.add_argument("--sliceSize", help="ADVANCED: number of bytes to process at a time in each thread", type=int, action = "store", default = 1000000)
    parser.add_argument("--verbose", help="Verbose output.", action = "store_true")
    parser.add_argument("--test", help="Test - runs 10 slices", action='store_true')


    args = parser.parse_args()

    ############## parse populations

    popNames = []
    popInds = []
    for p in args.population:
        popNames.append(p[0])
        if len(p) > 1: popInds.append(p[1].split(","))
        else: popInds.append([])

    if args.popsFile:
        with open(args.popsFile, "rt") as pf: popDict = dict([ln.split() for ln in pf])
        for ind in popDict.keys():
            try: popInds[popNames.index(popDict[ind])].append(ind)
            except: pass

    for p in popInds: assert len(p) >= 1, "All populations must be represented by at least one sample."

    allInds = list(set([i for p in popInds for i in p]))

    if args.ploidy is not None:
        ploidy = args.ploidy if len(args.ploidy) != 1 else args.ploidy*len(allInds)
        assert len(ploidy) == len(allInds), "Incorrect number of ploidy values supplied."
        ploidyDict = dict(zip(allInds,ploidy))
    elif args.ploidyFile is not None:
        with open(args.ploidyFile, "rt") as pf: ploidyDict = dict([[s[0],int(s[1])] for s in [l.split() for l in pf]])
    else: ploidyDict = dict(zip(allInds,[2]*len(allInds)))

    if args.haploid:
        for indName in args.haploid: ploidyDict[indName] = 1

    sampleData = genomics.SampleData(popNames = popNames, popInds = popInds, ploidyDict = ploidyDict)


    ############################################################################################################################################

    #open files

    if args.genoFile:
        if args.genoFile[-3:] == ".gz": genoFile = gzip.open(args.genoFile, "rt")
        else: genoFile = open(args.genoFile, "rt")
    else: genoFile = sys.stdin

    headerLine= genoFile.readline()

    if args.outFile:
        if args.outFile[-3:] == ".gz": outFile = gzip.open(args.outFile, "wt")
        else: outFile = open(args.outFile, "wt")
    else: outFile = sys.stdout

    outFile.write("scaffold\tposition\t")
    outFile.write("\t".join(popNames) + "\n")

    ##########################################################################################################################

    asCounts = args.asCounts if args.target else True
    keepNanLines = args.keepNanLines if args.target else True
    minData = args.minData if args.target else 0

    ##########################################################################################################################

    #counting stat that will let keep track of how far we are
    slicesQueued = 0
    resultsReceived = 0
    resultsWritten = 0
    linesWritten = 0

    '''Create queues to hold the data one will hold the line info to be passed to the analysis'''
    inQueue = SimpleQueue()
    #one will hold the results (in the order they come)
    resultQueue = SimpleQueue()
    #one will hold the sorted results to be written
    writeQueue = SimpleQueue()


    '''start worker Processes for analysis. The command should be tailored for the analysis wrapper function
    of course these will only start doing anything after we put data into the line queue
    the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
    workerThreads = []
    sys.stderr.write("\nStarting {} worker threads\n".format(args.threads))
    for x in range(args.threads):
        workerThread = Process(target=freqs_wrapper, args = (inQueue, resultQueue, headerLine, args.genoFormat, sampleData,
                                                        args.target, minData, asCounts, args.threshold, keepNanLines,))
        workerThread.daemon = True
        workerThread.start()
        workerThreads.append(workerThread)


    '''thread for sorting results'''
    sorterThread = Thread(target=sorter, args=(resultQueue,writeQueue,args.verbose, args.threads,))
    sorterThread.daemon = True
    sorterThread.start()

    '''start thread for writing the results'''
    writerThread = Thread(target=writer, args=(writeQueue, outFile, args.verbose,))
    writerThread.daemon = True
    writerThread.start()

    '''start background Thread that will run a loop to check run statistics and print
    We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
    checkerThread = Thread(target=checkStats)
    checkerThread.daemon = True
    checkerThread.start()

    ########################################################################################################################

    #generate slices and queue
    fileSlices = fileSlicer(genoFile, args.sliceSize)

    if not args.test:
        for fileSlice in fileSlices:
            inQueue.put((slicesQueued,fileSlice))
            slicesQueued += 1
    else:
        for fileSlice in fileSlices:
            inQueue.put((slicesQueued,fileSlice))
            slicesQueued += 1
            if slicesQueued == 10: break


    ############################################################################################################################################

    #Now we send completion signals to all worker threads
    for x in range(args.threads):
        inQueue.put((-1,None,)) # -1 tells the threads to break

    sys.stderr.write("\nClosing worker threads\n".format(args.threads))
    for x in range(len(workerThreads)):
        workerThreads[x].join()

    sorterThread.join()
    writerThread.join()

    sys.stderr.write("\nDone\n")

    genoFile.close()
    outFile.close()

    sys.exit()

