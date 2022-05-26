import gzip

import genomics


####################################################################################################################

genoFileName = "hel.18.DP10GQ30min10.BI.geno.gz"
genoFile = gzip.open(genoFileName, "r")

#samples

popNames = ["ama","agl","melG","thxn","ser"]

samples = [["ama.JM160","ama.JM216","ama.JM293","ama.JM48"],
           ["agl.JM108","agl.JM112","agl.JM569","agl.JM572"],
           ["melG.CJ13435","melG.CJ9315","melG.CJ9316","melG.CJ9317"],
           ["thxn.JM313","thxn.JM57","thxn.JM84","thxn.JM86"],
           ["ser.JM202"]]


sampleData = genomics.SampleData(popInds = samples, popNames = popNames)

P1 = "melG"
P2 = "ama"
P3 = "thxn"
O = "ser"


#####################################################################################################


outFileName = "ABBA_BABA_melG_ama_thxn_ser.w50m1s25.csv"
outFile = open(outFileName, "w")

outFile.write("scaffold,start,end,sites,ABBA,BABA,D,fd\n")

#####################################################################################################

#get windows and analyse
windSize = 50000
stepSize = 25000
minSites = 1000

windowGenerator = genomics.slidingWindows(genoFile, windSize, stepSize, )
n=0
for window in windowGenerator:
    if window.seqLen() >= minSites:
        #make alignment object
        Aln = genomics.callsToAlignment(window.seqDict(), sampleData, seqType = "pairs")
        #get ABBA BABA stats
        statsDict = genomics.ABBABABA(Aln, P1, P2, P3, O)
        stats = [statsDict[s] for s in ["ABBA","BABA","D","fd"]]
        #add window stats and write to file
        out_data = [window.scaffold, window.start, window.end, window.seqLen()] + [round(s,3) for s in stats]
        outFile.write(",".join([str(x) for x in out_data]) + "\n")
    print n
    n+=1


outFile.close()
genoFile.close()

