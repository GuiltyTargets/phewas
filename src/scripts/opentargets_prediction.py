# -*- coding: utf-8 -*-

"""Calculates the AUC for predictions of drug target proteins."""

# part 7a

from collections import defaultdict
import logging
import os

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from guiltytargets.ppi_network_annotation.parsers import parse_association_scores
from guiltytargets_phewas.utils import timed_main_run
from guiltytargets_phewas.constants import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('p7b.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

assert os.path.isdir(DATA_BASE_DIR), "Update your data_base_dir folder for this environment."

"""
hc
AUC-PR:   0.47087534303574397
AUC-ROC:  0.8950764506647462
lc
AUC-PR:   0.7062949094686672
AUC-ROC:  0.9735694536274568
aml
AUC-PR:   0.37372927153167096
AUC-ROC:  0.8908891379870963
ipf
AUC-PR:   0.46839457215896996
AUC-ROC:  0.933710525259821
ms
AUC-PR:   0.7915300099433009
AUC-ROC:  0.9633356773318599
ad
AUC-PR:   0.4445401652632658
AUC-ROC:  0.963198301702073
"""

# Paths
dge_base_path = os.path.join(DATA_BASE_DIR, 'DGE')
targets_base_path = os.path.join(DATA_BASE_DIR, 'OpenTargets')
string_base_path =os.path.join(DATA_BASE_DIR, 'STRING')

g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'linkprediction')
string_graph_path = os.path.join(string_base_path, 'string_entrez.edgelist')

max_log2_fold_change, min_log2_fold_change = -1 * lfc_cutoff, lfc_cutoff


# TODO Move these functions to a common file to be used by the scripts.
def targets_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_entrez.txt')


def assoc_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_assoc_entrez.txt')


def main():
    """ """
    results = defaultdict(list)
    for disease_abv in ['hc', 'lc', 'aml', 'ipf', 'ms', 'ad']:
        # disease_abv =
        assoc_path = assoc_file(disease_abv)
        assoc_scores = parse_association_scores(assoc_path)

        targets_path = targets_file(disease_abv)
        df = pd.read_csv(
            targets_path,
            names=['gene'],
            dtype={'gene': str},
            sep='\t'
        )
        genes = list(df['gene'])
        df = pd.DataFrame(
            [
                (x, y, 1 if x in genes else 0)
                for x, y
                in assoc_scores.items()
            ]
        )
        results['dge'].append(disease_abv)
        results['aps'].append(average_precision_score(df[2], df[1]))
        results['auc'].append(roc_auc_score(df[2], df[1]))
    pd.DataFrame(results).to_csv(os.path.join(g2v_path, 'ot_target_prediction.tsv'), sep='\t')


if __name__ == '__main__':
    timed_main_run(main, logger)
