# -*- coding: utf-8 -*-

"""Compares the weighting of samples using the association score for the  to unweighted approach. Uses gat2vec
default parameters and logistic regression."""

# part 3

import logging
import os
import warnings

# Suppress warnings
warnings.simplefilter('ignore')

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_association_scores, parse_gene_list
from guiltytargets_phewas.utils import timed_main_run
from guiltytargets_phewas.constants import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('p3_w_vs_unw.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

assert os.path.isdir(DATA_BASE_DIR), "Update your data_base_dir folder for this environment."

# Paths
dge_base_path = os.path.join(DATA_BASE_DIR, 'DGE')
targets_base_path = os.path.join(DATA_BASE_DIR, 'OpenTargets')
string_base_path =os.path.join(DATA_BASE_DIR, 'STRING')

g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'part3')
phewas_path = os.path.join(DATA_BASE_DIR, 'phewas_catalog', 'phewas_entrez.txt')
string_graph_path = os.path.join(string_base_path, 'string_entrez.edgelist')

min_log2_fold_change, max_log2_fold_change = -1 * lfc_cutoff, lfc_cutoff


# TODO Move these functions to a common file to be used by the scripts.
def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def targets_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_entrez.txt')


def assoc_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_assoc_entrez.txt')


def dge_file(dge_code: str) -> str:
    file = 'DifferentialExpression' + ('.csv' if dataset_to_disease_abv(dge_code) == 'ad' else '.tsv')
    return os.path.join(dge_base_path, dge_code, file)


disease_efo_dict = {
    'ad': 'EFO_0000249',
    'hc': 'EFO_0000182',
    'aml': 'EFO_0000222',
    'ipf': 'EFO_0000768',
    'lc': 'EFO_0001422',
    'ms': 'EFO_0003885'
}


def get_ppi_results(
    ppi_graph_path: str,
        dataset: str,
        evaluation: str = 'cv',
        assoc_path: str = None
) -> pd.DataFrame:
    dge_params = dge_params_ad if dataset in AD_DGE_DATASETS else dge_params_dis
    gene_list = parse_dge(
        dge_path=dge_file(dataset),
        entrez_id_header=dge_params['id'],
        log2_fold_change_header=dge_params['l2f'],
        adj_p_header=dge_params['adjp'],
        entrez_delimiter=split_char,
        base_mean_header=dge_params['mean'],
    )
    network = generate_ppi_network(
        ppi_graph_path=ppi_graph_path,
        dge_list=gene_list,
        max_adj_p=max_padj,
        max_log2_fold_change=max_log2_fold_change,
        min_log2_fold_change=min_log2_fold_change,
        ppi_edge_min_confidence=ppi_edge_min_confidence,
        current_disease_ids_path='',
        disease_associations_path=phewas_path,
    )

    logger.info(f'Nodes {len(network.graph.vs)}')
    targets = parse_gene_list(targets_file(dataset_to_disease_abv(dataset)), network)
    assoc_score = assoc_path and parse_association_scores(assoc_path)
    logger.debug(f'Number of targets being used for the network: {len(targets)}')

    write_gat2vec_input_files(
        network=network,
        targets=targets,
        home_dir=g2v_path,
        assoc_score=assoc_score
    )
    metrics_df, _ = rank_targets(
        directory=g2v_path,
        network=network,
        evaluation=evaluation,
    )
    df = pd.DataFrame()
    df['auc'] = metrics_df['auc']
    df['aps'] = metrics_df['aps']
    df['eval'] = 'weighted' if assoc_path else 'unweighted'
    df['dge'] = dataset
    return df


def main():
    """ """

    for key, dataset in DGE_DATASETS.items():
        results = pd.DataFrame()
        for ds in dataset:
            # Weighted
            part_df = get_ppi_results(
                string_graph_path,
                ds,
                evaluation='nested_cv',
                assoc_path=assoc_file(dataset_to_disease_abv(ds))
            )
            results = results.append(part_df, ignore_index=True)
            # Unweighted
            part_df = get_ppi_results(
                string_graph_path,
                ds,
                evaluation='nested_cv'
            )
            results = results.append(part_df, ignore_index=True)
        results.to_csv(os.path.join(g2v_path, f'results_{key}_df.tsv'), sep='\t')
        # Prepare results
        for metric in ['auc', 'aps']:
            fig = sns.boxplot(
                data=results,
                x='dge',
                y=metric,
                hue='eval'
            ).get_figure()
            fig.suptitle('Comparison of weighting negative samples vs unweighted.')
            fig.savefig(f'comp3_{key}_{metric}.png')
            plt.close()


if __name__ == '__main__':
    timed_main_run(main, logger)
