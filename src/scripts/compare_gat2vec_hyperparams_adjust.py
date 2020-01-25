# -*- coding: utf-8 -*-

"""Analyzes the gat2vec parameters to find optimal non-default values."""

# part 4

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
fh = logging.FileHandler('p4_g2v_hyperparam_opt.log')
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

g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'part4')
phewas_path = os.path.join(DATA_BASE_DIR, 'phewas_catalog', 'phewas_entrez.txt')
string_graph_path = os.path.join(string_base_path, 'string_entrez.edgelist')

max_log2_fold_change, min_log2_fold_change = -1 * lfc_cutoff, lfc_cutoff


# TODO Move these functions to a common file to be used by the scripts.
def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def targets_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_entrez.txt')


def assoc_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_assoc_entrez.txt')


def dge_file(dge_code: str) -> str:
    ext = '.csv' if dataset_to_disease_abv(dge_code) == 'ad' else '.tsv'
    return os.path.join(dge_base_path, dge_code, 'DifferentialExpression' + ext)


disease_efo_dict = {
    'ad': 'EFO_0000249',
    'hc': 'EFO_0000182',
    'aml': 'EFO_0000222',
    'ipf': 'EFO_0000768',
    'lc': 'EFO_0001422',
    'ms': 'EFO_0003885'
}


def optimize_g2v_parameters(
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
        csv_separator=';'
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
    results = pd.DataFrame()
    for dim in [32, 64, 128, 256]:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=evaluation,
            dimension=dim,
        )
        df = pd.DataFrame()
        df['auc'] = metrics_df['auc']
        df['aps'] = metrics_df['aps']
        df['eval'] = str(dim)
        df['dge'] = dataset
        df['param'] = 'Dimension'
        logger.debug('df')
        logger.debug(df)
        results = results.append(
            assemble_results_df(metrics_df, dataset, 'Dimension', dim),
            ignore_index=True
        )
    """num_walks: int = 10
    walk_length: int = 80
    dimension: int = 128
    window_size: int = 5"""
    for nw in [6, 10, 20, 40, 80]:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=evaluation,
            dimension=256,
            num_walks=nw,
        )
        results = results.append(
            assemble_results_df(metrics_df, dataset, 'Num Walks', nw),
            ignore_index=True
        )
    for wl in [20, 40, 80, 120, 160]:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=evaluation,
            dimension=256,
            walk_length=wl
        )
        results = results.append(
            assemble_results_df(metrics_df, dataset, 'Walk Length', wl),
            ignore_index=True
        )
    for ws in [3, 5, 7, 10, 20, 40]:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=evaluation,
            dimension=256,
            window_size=ws
        )
        results = results.append(
            assemble_results_df(metrics_df, dataset, 'Window Size', ws),
            ignore_index=True
        )
    logger.debug('results')
    logger.debug(results)
    return results


def assemble_results_df(
        metrics_df: pd.DataFrame,
        dataset: str,
        param_name: str,
        param_val: int
) -> pd.DataFrame:
    """Creates a partial dataframe to be stored for plotting.

    :param metrics_df: A dataframe with the important metrics.
    :param dataset: The name of the dataset.
    :param param_name: The hyperparameter being optimized in the metrics dataframe.
    :param param_val: The value of the hyperparameter being optimized.
    :return: A dataframe for plotting.
    """
    df = pd.DataFrame()
    df['auc'] = metrics_df['auc']
    df['aps'] = metrics_df['aps']
    df['eval'] = str(param_val)
    df['dge'] = dataset
    df['param'] = param_name
    return df


def main():
    """ """
    results = pd.DataFrame()
    for dataset in DGE_DATASETS.values():
        for ds in dataset:
            # Weighted
            part_df = optimize_g2v_parameters(
                string_graph_path,
                ds,
                evaluation='nested_cv',
            )
            results = results.append(part_df, ignore_index=True)
    results.to_csv(os.path.join(g2v_path, f'results_df.tsv'), sep='\t')


if __name__ == '__main__':
    timed_main_run(main, logger)