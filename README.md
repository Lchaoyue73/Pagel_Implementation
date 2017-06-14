# Manual

This is a manual for the implementation of calculating the Pagel likelihood ratios and running hierarchical clustering.


## Prerequisites

  - [Python 2.7](https://www.python.org/downloads/)

  - Pagel's Software: [Bayestraits](http://www.evolution.rdg.ac.uk/BayesTraitsV3/BayesTraitsV3.html). Download the appropriate version for your system, and move the executable application 'BayestraitsV3' to the same folder as python scripts.

  - Python libraries: 
  	  - [DendroPy Phylogenetic Computing Library](https://www.dendropy.org)
      - Numpy, Pandas, Scipy

  - A csv file containing phylogenetic profiles (see below)
  - A phylogenetic tree file (see below). 

#### Phylogenetic profile matrix

A csv file ("," seperated) containing phylogenetic profiles is required. Each row is a profile for one gene across given genomes.

|       | genome1 | genome2 | genome3 | genome4 | ... |
|:-----:|:-------:|:-------:|:-------:|:-------:|:---:|
| gene1 |    0    |    1    |    0    |    0    | ... |
| gene2 |    0    |    0    |    0    |    1    | ... |
| gene3 |    1    |    0    |    0    |    1    | ... |
| gene4 |    1    |    1    |    0    |    1    | ... |
| gene5 |    0    |    1    |    0    |    1    | ... |
|  ...  |   ...   |   ...   |   ...   |   ...   | ... |

#### Phylogenetic tree

A phylogenetic tree file (nexus format) is required. The tree format needs to satify the requirements of Pagel's Bayestraits software, which is that the taxa names must not be included in the descriptionof the tree but should be linked to a number in the translate section of the tree file. Details can be found in the [Bayestraits mannual](http://www.evolution.rdg.ac.uk/BayesTraitsV3/Files/BayesTraitsV3.Manual.pdf). 

A [BayestreeConverter](http://www.evolution.rdg.ac.uk/BayesTrees.html), which is also provided by Pagel, can be used to generate the right tree format.


#### Command files

The command files are used to contain the commands of Bayestraits to run. Two default command files are included in the directory. commandfile0.txt is for independent evolution model and commandfile1.txt is for dependent evolution model.Details can be found in the [Bayestrait mannual](http://www.evolution.rdg.ac.uk/BayesTraitsV3/Files/BayesTraitsV3.Manual.pdf). The command files must be located in the same directory as python scripts. 



## Example 

#### Step 1: Calculate the likelihood ratio statistics

- arguments: 
    - -h : help
	- -t : a phylogenetic tree
	- -m : phylogenetic profile matrix

- e.g.

	python CalculateStatistics.py -t test_tree.trees -m test_profile.csv


- output: LR_outputs.txt with columns: gene1, gene2, indepdent, dependent and likelihood ratio.


#### Step 2: Run hierarchical clustering

The LR_outputs.txt must be under the same directory.

- arguments:
    - h : help
	- v : a numeric value - tHe therreshold to apply when forming clusterings.
	- c : the criterion to use in forming clusters. 
	See Details in the document for [scipy.cluster.hierarchy.fcluster](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster).

- e.g.

	python RunClustering.py -v 7 -c distance


- output: Clustering_outputs.txt. Two columns in this file: the first column is the list of input genes; the second column is the list of corresponding cluster labels.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

