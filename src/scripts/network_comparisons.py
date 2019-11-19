# -*- coding: utf-8 -*-

"""Compares the GuiltyTargets approach using different input networks: STRING, HumanBase and OpenBEL. Uses gat2vec
default parameters and logistic regression."""

# part 1

import logging
import os
import warnings

# Suppress warnings and BEL parsing errors
warnings.simplefilter('ignore')
# Suppress Pybel parsing errors
logging.getLogger('pybel.parser').setLevel(logging.CRITICAL)
logging.getLogger('pybel_tools.input.reified_graph.input').setLevel(logging.CRITICAL)
debug = logging.getLogger('pybel.parser')

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_gene_list
from guiltytargets_phewas.constants import *
from guiltytargets_phewas.utils import generate_bel_network, timed_main_run

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('p1_network_comparisons.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

assert os.path.isdir(DATA_BASE_DIR), "Update your data_base_dir folder for this environment."

# Paths
dge_base_path = os.path.join(DATA_BASE_DIR, 'DGE')
tissue_base_path = os.path.join(DATA_BASE_DIR, 'HumanBase')
targets_base_path = os.path.join(DATA_BASE_DIR, 'OpenTargets')
string_base_path =os.path.join(DATA_BASE_DIR, 'STRING')

g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'part1')
bel_files = os.path.join(DATA_BASE_DIR, 'bel_data', 'all')  # 126 targets
string_graph_path = os.path.join(string_base_path, 'string_entrez.edgelist')

max_log2_fold_change, min_log2_fold_change = -1 * lfc_cutoff, lfc_cutoff

network_alias = {
    'cerebral_cortex_top': 'ctex',
    'hippocampus_top': 'hipp',
    'liver_top': 'hb',
    'bone_marrow_top': 'hb',
    'lung_top': 'hb',
    'central_nervous_system_top': 'hb',
    'string_entrez.edgelist': 'st'
}

# igraph.read requires the files to be gunzipped. Have to do that after downloading.
hb_ppi = {
    'BM10': ['cerebral_cortex_top', 'hippocampus_top'],
    'BM22': ['cerebral_cortex_top', 'hippocampus_top'],
    'BM36': ['cerebral_cortex_top', 'hippocampus_top'],
    'BM44': ['cerebral_cortex_top', 'hippocampus_top'],
    'Mayo_CTX': ['cerebral_cortex_top', 'hippocampus_top'],
    'Mayo_CBE': ['cerebral_cortex_top', 'hippocampus_top'],
    'RosMap': ['cerebral_cortex_top', 'hippocampus_top'],
    'hc': ['liver_top'],
    'aml': ['bone_marrow_top'],
    'ipf': ['lung_top'],
    'lc': ['liver_top'],
    'ms': ['central_nervous_system_top']
}


# TODO Move these functions to a common file to be used by the scripts.
def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def dge_file(dge_code: str) -> str:
    ext = '.csv' if dataset_to_disease_abv(dge_code) == 'ad' else '.tsv'
    return os.path.join(dge_base_path, dge_code, 'DifferentialExpression' + ext)


def targets_file(disease):
    return os.path.join(targets_base_path, disease, 'ot_entrez_missing.txt')


def get_bel_results(dataset) -> (pd.DataFrame, pd.DataFrame):
    """"""
    dge_params = dge_params_ad if dataset in AD_DGE_DATASETS else dge_params_dis
    gene_list = parse_dge(
        dge_path=dge_file(dataset),
        entrez_id_header=dge_params['id'],
        log2_fold_change_header=dge_params['l2f'],
        adj_p_header=dge_params['adjp'],
        entrez_delimiter=split_char,
        base_mean_header=dge_params['mean'],
    )
    network = generate_bel_network(
        bel_graph_path=bel_files,
        dge_list=gene_list,
        max_adj_p=max_padj,
        max_log2_fold_change=max_log2_fold_change,
        min_log2_fold_change=min_log2_fold_change,
    )
    logger.info(f'Nodes {len(network.graph.nodes)}')
    targets = parse_gene_list(
        os.path.join(targets_base_path, dataset_to_disease_abv(dataset), 'ot_symbol.txt'),
        network
    )
    logger.debug(f'Number of targets being used for the network: {len(targets)}')

    write_gat2vec_input_files(
        network=network,
        targets=targets,
        home_dir=g2v_path,
    )
    metrics_df, _ = rank_targets(
        directory=g2v_path,
        network=network,
        num_walks=30,
        walk_length=4,
        dimension=256,
        window_size=10
    )
    df = pd.DataFrame()
    df['auc'] = metrics_df['auc']
    df['aps'] = metrics_df['aps']
    df['ds'] = 'bel'
    df['dge'] = dataset
    return df


# TODO Refactor the get_ppi_results from all scripts.
def get_ppi_results(ppi_graph_path: str, dataset: str) -> pd.DataFrame:
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
    )

    logger.info(f'Nodes {len(network.graph.vs)}')
    logger.info(f'Edges {len(network.graph.es)}')
    targets = parse_gene_list(targets_file(dataset_to_disease_abv(dataset)), network)
    logger.debug(f'Number of targets being used for the network: {len(targets)}')

    write_gat2vec_input_files(
        network=network,
        targets=targets,
        home_dir=g2v_path,
    )
    metrics_df, _ = rank_targets(
        directory=g2v_path,
        network=network
    )
    df = pd.DataFrame()
    df['auc'] = metrics_df['auc']
    df['aps'] = metrics_df['aps']
    df['ds'] = network_alias[os.path.basename(ppi_graph_path)]
    df['dge'] = dataset
    return df


def main():
    """ """
    results = pd.DataFrame()
    for datasets in DGE_DATASETS.values():
        for ds in datasets:
            # OpenBEL
            if ds in AD_DGE_DATASETS:
                part_df = get_bel_results(ds)
                results = results.append(part_df, ignore_index=True)

            # STRING
            part_df = get_ppi_results(string_graph_path, ds)
            results = results.append(part_df, ignore_index=True)

            # HumanBase
            for tissue_file in hb_ppi[ds]:
                part_df = get_ppi_results(os.path.join(tissue_base_path, tissue_file), ds)
                results = results.append(part_df, ignore_index=True)
    results.to_csv(os.path.join(g2v_path, f'results_df.tsv'), sep='\t')

    # Prepare results
    for metric in ['auc', 'aps']:
        for key, datasets in DGE_DATASETS.items():
            fig = sns.boxplot(
                data=results[results['dge'].isin(datasets)],
                x='dge',
                y=metric,
                hue='ds'
            ).get_figure()
            fig.suptitle("Network comparison: STRING, BEL, HumanBase")
            fig.savefig(f'comp1_{metric}_{key}.png')  # AUC non-AD and AD
            plt.close()


if __name__ == '__main__':
    timed_main_run(main, logger)
