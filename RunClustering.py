import pandas as pd
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform
import numpy as np
import argparse


# Transform output file to a distance matrix
def LongToMatrix(geneList, lrFile):
    disDF = pd.DataFrame(0, index=geneList, columns=geneList)
    for i in range(len(geneList)):
        gi = geneList[i]
        giRows = lrFile[lrFile.gene1 == gi][['gene2', 'lr']]
        if giRows.shape[0] > 0:
            disDF.loc[gi, giRows.gene2] = giRows.lr.values
    return disDF


# likelihood ratio matrix to distance matrix
def ToDistanceMatrix(disDF):
    disMatrix = disDF.as_matrix()
    disMatrix = disMatrix + np.transpose(disMatrix)
    disMatrix = disMatrix.max() - disMatrix
    np.fill_diagonal(disMatrix, 0)
    return disMatrix


# Construct the condensed distance matrix for linkage function
def ToCondensedMatrix(disMatrix):
    # condense data format for hierarchical clustering
    disArray = squareform(disMatrix)
    return disArray


# Run hierarchical clustering
def RunHclust(t, criterior):
    z = hcl.linkage(disArray, method='average')
    clusterLabels = hcl.fcluster(z, t, criterion='distance')
    return clusterLabels


# Parse the inputs from command line
parser = argparse.ArgumentParser(description='Run Hierarchical Clustering')
parser.add_argument("-c", dest="criterion", type=str,
                    help="choose one of criteria: distance, maxclust, inconsistent...")
parser.add_argument("-v", dest="value", type=float,
                    help="value for criterion")
args = parser.parse_args()
t = args.value
criterion = args.criterion

# Read Likelihood ratio output file
lrFile = pd.read_csv("./LR_outputs.txt", names=["gene1", "gene2", "indep", "dep", "lr"])

# Prepare the data for hierarchical clustering
geneList = np.append(lrFile.gene2.unique(), lrFile.gene1.iloc[-1])
disDF = LongToMatrix(geneList, lrFile)
disMatrix = ToDistanceMatrix(disDF)
disArray = ToCondensedMatrix(disMatrix)

# Run clustering and output
clusterLabels = RunHclust(t, criterion)
clusterResults = pd.DataFrame({'GI': geneList, 'Labels': clusterLabels})
clusterResults.to_csv("./Clustering_outputs.txt", index=False)
