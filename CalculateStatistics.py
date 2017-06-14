import argparse
import os
import dendropy
import pandas as pd
import subprocess


# Check if files are valid
def is_valid_file(parser, inputFile):
    if not os.path.exists(inputFile):
        parser.error("The file %s does not exist." % inputFile)
    else:
        return open(inputFile, 'r')


# Extract the profile of a given pair of genes
def GetPairProfiles(gene1, gene2, fullProfile, tree, write=False):
    # get all taxa on the tree
    subGenomes = [tx.label for tx in tree.taxon_set]
    # get the profile and filenames
    pairProfile = fullProfile.ix[[gene1, gene2], subGenomes]
    pairProfile = pairProfile.transpose()
    filename = 'pairProfile_' + str(gene1) + '_' + str(gene2) + '_' + '.txt'
    if write:
        pairProfile.to_csv(filename, sep='\t', header=False, index=True)
    return filename


# Calculate the likelihood for a pair by calling Bayestraits
def RunPagelPair(treeName, pairfilename, commandN):
    commandfile = 'commandfile' + str(commandN) + '.txt'
    logfile = pairfilename + '.log.txt'
    subprocess.call(["./BayesTraitsV2", treeName, pairfilename],
                    stdin=open(commandfile, 'r'), shell=False)
    with open(logfile, 'r') as f:
            for line in f:
                pass
    lh = line.split()[1]
    os.remove(logfile)
    return lh


# Write results into output file
def WriteToOutputs(resultsfile, gene1, gene2, indep, dep, lr):
    with open(resultsfile, 'a') as results:
        results.write(",".join(map(str, [gene1, gene2, indep, dep, lr])) + "\n")


parser = argparse.ArgumentParser(description='Calculate Pagel likelhood')
parser.add_argument('-t', dest="tree", required=True,
                    help="Phylogenetic tree file (nexus format)",
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument('-m', dest="profile", required=True,
                    help="Phylogentic profile matrix",
                    type=lambda x: is_valid_file(parser, x))
args = parser.parse_args()
fullProfile = pd.read_csv(args.profile, index_col=0)
tree = dendropy.Tree.get_from_stream(args.tree, 'nexus')
treeName = args.tree.name
genomes = [tx.label for tx in tree.taxon_set]
resultsfile = 'LR_outputs.txt'
# pdb.set_trace()
geneList = fullProfile.index.values


# Remove the old output file
try:
    os.remove(resultsfile)
except OSError:
    pass


# Run all pairs in the gene list
for ith in range(len(geneList)):
    gene1 = geneList[ith]
    for jth in range(0, ith):
        gene2 = geneList[jth]
        pairfilename = GetPairProfiles(gene1, gene2, fullProfile, tree, write=True)
        # Run the indepent (0) and dependent (1) models
        indep = RunPagelPair(treeName, pairfilename, 0)
        dep = RunPagelPair(treeName, pairfilename, 1)
        # Calculate the likelihood
        lr = 2 * round(float(dep) - float(indep), 6)
        # Some pairs might have small negtive likelihood ratios
        if lr < 0:
            lr = 0
        WriteToOutputs(resultsfile, gene1, gene2, indep, dep, lr)
        os.remove(pairfilename)
