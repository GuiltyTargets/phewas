GuiltyTargets Extension
=======================
This is a tool for therapeutic target prioritization using network representation learning. While the original [GuiltyTargets](https://github.com/GuiltyTargets/guiltytargets) focused on established PPI networks like STRING and HIPPIE, this work explores possibilities like OpenBEL documents and Tissue Specific networks as well as disease-gene relationships.

Under development: Link Prediction for target repositioning.

Installation
------------
Download this repository, go to the directory it resides and run:

.. code-block:: bash

   $ git clone https://github.com/guiltytargets/guiltytargets_phewas.git
   $ cd guiltytargets_phewas
   $ pip install -e .

Usage
-----
After that, you can use the python scripts (found in src/scripts) to check the obtained results:

1. network_comparisons.py

In this script, different Protein networks for GuiltyTargets were compared to each other: the original approach,
using the STRING PPI network; a reified OpenBEL Graph and tissue specific PPI networks from the GIANT webservice.

2. compare_svm_rr.py

The use of Biased SVM is compared against penalized Ridge Regression. Those
methods run under nested cross validation for the proper optimization of
class weights and the penalization hyperparameters.

3. compare_weighted_unweighted.py

The association score between genes and the analyzed disease is used for
sample weighting. In this case, the positive targets are given maximum sample
weights (one) and the negative samples are being weighted with
1 - association_score, giving less certainty about the outcome for highly
associated targets.

4. compare_gat2vec_hyperparams_adjust.py

The hyper parameters of gat2vec (walk length, window size, number of walks
and dimension) are optimized for the network.

5. link_prediction.py

A different approach to drug reposition is analyzed: building an heterogeneous
graph containing associations among drug, disease and targets and trying to
predict missing links between the disease and its targets.

INPUT FILES (WIP)
-----------
There are 5 files which are necessary to run this program. All input files should be found under input_directory 

1. ``ppi_graph_path``:  A path to a file containing a protein-protein interaction network in the format of:

    **EntrezID** **EntrezID** **CONFIDENCE**
    
    
    Such as:
    
    216 216 0.76
    
    3679 1134 0.73
    
    55607 71 0.65
    
    5552 960 0.63
    
    2886 2064 0.9
    
    5058 2064 0.73
    
    1742 2064 0.87
    
    An example of such a network can be found [here](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)

2. ``dge_path``: A path to a file containing an experiment, in tsv format. Rows show individual entries, columns are the values of the following properties:
  - **Base mean**
  - **Log fold change**
  - **Adjusted p value**
  - **Entrez id**

  The file may contain other columns too, but the indices and names of the above columns must be entered to the configuration file.

3. ``targets_path``: A path to a file containing a list of Entrez ids of known targets, in the format of

    EntrezID1
    
    EntrezID2
    
    ...
    
    
    Such as:
    
    1742
    
    3996
    
    150
    
    152
    
    151

4. ``assoc_path``: A path to a file containing a list of Entrez ids and the known association scores to the current investigated disease, in the format of:

    EntrezID1	SCORE1
    
    EntrezID2	SCORE2
    
    ...
    
    
    Such as:
    
    1742	0.7
    
    3996	0.2
    
    2150	0.1
    
    3152	0.0
    
    5151	0.8

5. ``phewas_path``: A path to a file containing a list of Entrez ids and the known association scores to the current investigated disease, in the format of:

OPTIONS
-------
The options that should be set are:

max_adj_p: Maximum value for adjusted p-value for a gene to be considered differentially expressed.

entrez_delimiter: If there is more than one Entrez id per row in the diff. expr. file, the separator betweem them.
