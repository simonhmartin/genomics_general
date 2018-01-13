
#code to parse multiple trees from ms or msms
# outputs a cleaned tree file (without the "[]" giving number of sites) and a tree positions tsv file.

import sys, argparse
import numpy as np
import gzip


def merge_identical_topos(trees):
    output = {"trees":[], "indices":[]}
    n=0
    last = None
    for tree in trees:
        #for each tree, get the unique topo ID
        tree.unroot()
        ID = tree.get_topology_id()
        if ID!=last:
            #its unique
            output["trees"].append(tree)
            output["indices"].append([n])
            last=ID
        else:
            #it matches its sister, so no need to add it, just get it's index
            output["indices"][-1].append(n)
        n+=1
    return output


######################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--msOutput", help="output from ms", action = "store", required = False)
parser.add_argument("-t", "--outTrees", help="output trees file", action = "store", required = False)
parser.add_argument("-d", "--outData", help="output data file", action = "store", required = False)
parser.add_argument("--remove_branch_lengths", action = "store_true")
parser.add_argument("--merge_identical_topologies", action = "store_true")


args = parser.parse_args()

if args.msOutput:
    if args.msOutput[-3:] == ".gz": msOutput = gzip.open(args.msOutput, "r")
    else: msOutput = open(args.msOutput, "r")
else: msOutput = sys.stdin

if args.outTrees:
    if args.outTrees[-3:] == ".gz": outTrees = gzip.open(args.outTrees, "w")
    else: outTrees = open(args.outTrees, "w")
else: outTrees = sys.stdout

if args.outData:
    if args.outData[-3:] == ".gz": outData = gzip.open(args.outData, "w")
    else: outData = open(args.outData, "w")

rmLens = args.remove_branch_lengths

merge = args.merge_identical_topologies
if merge: rmLens = True

######################################################################################

#get trees from ms output
msTrees = [l.strip() for l in msOutput if l.strip()[-1]==";"]

#get clean trees
trees = [t[t.index("("):] for t in msTrees]

if merge or rmLens:
    from ete3 import Tree

    trees = [Tree(tree) for tree in trees]

    if merge:
        #merge identical adjacent topologies
        merged = merge_identical_topos(trees)
        trees = merged["trees"]
        

    #text for trees
    if rmLens: trees_text = [tree.write(format=9) for tree in trees]
    else: trees_text = [tree.write(format=5) for tree in trees]

else: trees_text = trees

outTrees.write("\n".join(trees_text) + "\n")

outTrees.close()

######################################################################################

if args.outData:
    
    assert msTrees[0][0] == "[", "Site data not present."
    
    #get tree sites
    treeSites = np.array([int(t[:(t.index("]")+1)].strip("[]")) for t in msTrees])
    
    if merge:
        #if we've merged trees, we have to sum the sets of sites
        treeSites = np.array([np.sum(treeSites[indices]) for indices in merged["indices"]])
    
    #get start and end positions
    ends = np.cumsum(treeSites)
    starts = ends - treeSites + 1
    
    #write output
    outData.write("start\tend\n")
    for x in range(len(ends)): outData.write(str(starts[x]) + "\t" + str(ends[x]) + "\n")
    
    outData.close()

sys.exit()


