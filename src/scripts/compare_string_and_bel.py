# -*- coding: utf-8 -*-

"""Runs the ranking of drug targets using the HumanBase data as PPI network."""

import logging
import numpy as np
import os
import seaborn as sns
import time
import traceback
import warnings

import matplotlib.pyplot as plt
import pandas as pd
from pybel import from_json_file

from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_association_scores, parse_gene_list
from guiltytargets_phewas.utils import generate_bel_network

# from constants
DATA_BASE_DIR = r'c:/users/mauricio/thesis/data'

ppi_edge_min_confidence = 0.0
lfc_cutoff = 1.5

# for differential expression (AD files)
base_mean_name_ad = 'baseMean'
log_fold_change_name_ad = 'log2FoldChange'
adjusted_p_value_name_ad = 'padj'
entrez_id_name_ad = 'entrez'

# for differential expression (other files)
base_mean_name_dis = None
log_fold_change_name_dis = 'logFC'
adjusted_p_value_name_dis = 'adj.P.Val'
entrez_id_name_dis = 'Gene.ID'

# for differential expression (any)
max_padj = 0.05
split_char = '///'
diff_type = 'all'

# Suppress warnings and BEL parsing errors
warnings.simplefilter('ignore')
# Suppress Pybel parsing errors
logging.getLogger('pybel.parser').setLevel(logging.CRITICAL)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('optimize_parameters_ad.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

assert os.path.isdir(DATA_BASE_DIR), "Update your data_basedir folder for this environment."

dge = 'BM10'

# Paths
targets_file = os.path.join(DATA_BASE_DIR, 'OpenTargets', 'ad', 'ot_symbol.txt')
assoc_file = os.path.join(DATA_BASE_DIR, 'OpenTargets', 'ad', 'ot_assoc_symbol.txt')
g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'ppi', 'ad')
phewas_path = os.path.join(DATA_BASE_DIR, 'phewas_catalog', 'phewas_symbol.txt')

dge_path = os.path.join(DATA_BASE_DIR, 'DGE', 'AMP-AD', f'DifferentialExpression_{dge}.csv')
ppi_graph_path = os.path.join(DATA_BASE_DIR, 'STRING', 'string_symbol.edgelist')
bel_files = os.path.join(DATA_BASE_DIR, 'bel_data', 'all')  # 126 targets

min_log2_fold_change, max_log2_fold_change = -1 * lfc_cutoff, lfc_cutoff

assert os.path.isfile(targets_file)
assert os.path.isdir(g2v_path)
assert os.path.isfile(assoc_file)
assert os.path.isfile(phewas_path)

def do_rankings(
        network,
        targets,
        num_walks=10,
        walk_length=80,
        dimension=128,
        window_size=5,
):
    """

    :param network:
    :param targets:
    :param num_walks: For the ranking, uses the default value from gat2vec
    :param walk_length: For the ranking, uses the default value from gat2vec
    :param dimension: For the ranking, uses the default value from gat2vec
    :param window_size: For the ranking, uses the default value from gat2vec
    :return:
    """
    auc_df, aps_df = pd.DataFrame(), pd.DataFrame()
    targets = targets or parse_gene_list(targets_file, network)
    logger.debug(f'Number of targets being used for the STRING network: {len(targets)}')
    assoc_score = assoc_file and parse_association_scores(assoc_file)
    # Unweighted
    write_gat2vec_input_files(
        network=network,
        targets=targets,
        home_dir=g2v_path,
    )
    # TODO Use nested
    # for ev_method in ['cv', 'nested_cv', 'nested_svm']:
    for ev_method in ['cv', 'svm']:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=ev_method,
            num_walks=num_walks,
            walk_length=walk_length,
            dimension=dimension,
            window_size=window_size,
        )
        auc_df[ev_method] = metrics_df['auc']
        aps_df[ev_method] = metrics_df['aps']
    # Weighted
    write_gat2vec_input_files(
        network=network,
        targets=targets,
        home_dir=g2v_path,
        assoc_score=assoc_score,
    )
    # for ev_method in ['nested_cv', 'nested_svm']:
    for ev_method in ['cv', 'svm']:
        metrics_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            evaluation=ev_method,
            num_walks=num_walks,
            walk_length=walk_length,
            dimension=dimension,
            window_size=window_size,
        )
        label = 'w_' + ev_method
        auc_df[label] = metrics_df['auc']
        aps_df[label] = metrics_df['aps']
    return auc_df, aps_df


def prepare_plot(df_bel: pd.DataFrame, df_ppi: pd.DataFrame, title: str, out_file: str):
    df_bel['ds'] = 'OpenBEL'
    df_ppi['ds'] = 'PPI'
    df_join = df_ppi.append(df_bel, ignore_index=True)
    melt = df_join.melt(
        id_vars='ds',
        var_name='columns',
        # TODO Use nested
        value_vars=[
            'cv',
            'svm',
            'w_cv',
            'w_svm',
        ]
    )
    fig = sns.boxplot(y='value', x='columns', data=melt, hue='ds').get_figure()
    fig.suptitle(title)
    fig.savefig(out_file)


def get_string_results(targets) -> (pd.DataFrame, pd.DataFrame):
    gene_list = parse_dge(
        dge_path=dge_path,
        entrez_id_header=entrez_id_name_ad,
        log2_fold_change_header=log_fold_change_name_ad,
        adj_p_header=adjusted_p_value_name_ad,
        entrez_delimiter=split_char,
        base_mean_header=base_mean_name_ad,
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

    return do_rankings(network, targets)


def get_bel_results(targets) -> (pd.DataFrame, pd.DataFrame):
    gene_list = parse_dge(
        dge_path=dge_path,
        entrez_id_header=entrez_id_name_ad,
        log2_fold_change_header=log_fold_change_name_ad,
        adj_p_header=adjusted_p_value_name_ad,
        entrez_delimiter=split_char,
        base_mean_header=base_mean_name_ad,
    )
    network = generate_bel_network(
        bel_graph_path=bel_files,
        dge_list=gene_list,
        max_adj_p=max_padj,
        max_log2_fold_change=max_log2_fold_change,
        min_log2_fold_change=min_log2_fold_change,
        disease_associations_path=phewas_path
    )

    return do_rankings(
        network,
        targets,
        num_walks=30,
        walk_length=4,
        dimension=256,
        window_size=10,
    )


def main():
    targets_list = []
    auc_bel_df, aps_bel_df = get_bel_results(targets_list)
    auc_ppi_df, aps_ppi_df = get_string_results(targets_list)

    prepare_plot(df_bel=auc_bel_df, df_ppi=auc_ppi_df, title='AUROC results', out_file='comparison_str_bel_auc.png')
    prepare_plot(df_bel=aps_bel_df, df_ppi=aps_ppi_df, title='AUPRC results', out_file='comparison_str_bel_aps.png')


if __name__ == '__main__':
    start_time = time.time()
    logger.info('Starting...')
    try:
        main()
    except Exception as e:
        logger.error(type(e))
        logger.error(traceback.format_exc())
    finally:
        logger.info(f"Total time: {time.time() - start_time}")
