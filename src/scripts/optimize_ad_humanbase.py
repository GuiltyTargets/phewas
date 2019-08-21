# -*- coding: utf-8 -*-

"""Runs the ranking of drug targets using the HumanBase data as PPI network."""

import logging
import itertools
from math import ceil
import time
import traceback
import warnings

import matplotlib.pyplot as plt
import pandas as pd

from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_gene_list
from .constants import *

# Suppress warnings
warnings.simplefilter('ignore')

# Log the run
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
targets_file = os.path.join(DATA_BASE_DIR, r'OpenTargets/ad/ot_entrez.txt')
assoc_file = os.path.join(DATA_BASE_DIR, r'OpenTargets/ad/ot_assoc_entrez.txt')
g2v_path = os.path.join(DATA_BASE_DIR, r'gat2vec_files/ppi/ad')
phewas_path = None  # Phewas file need to be converted to Entrez

ppi_base_path = os.path.join(DATA_BASE_DIR, r'HumanBase')
dge_path = os.path.join(DATA_BASE_DIR, f'DGE/AMP-AD/DifferentialExpression_{dge}.csv')
graph_paths = ['hippocampus_top', 'cerebral_cortex_top']

def optimize_parameters(
        network,
        nwal: list,
        wlen: list,
        dim: list,
        ws_mult: list,
        evaluation: str='auc'
) -> int:
    """Tests the algorithm with the given parameters. Only one of the parameters can be tested at a time.

    :param network: Graph object.
    :param nwal: List of values to be used for number of walks.
    :param wlen: List of values to be used for walk length.
    :param dim: List of values to be used for dimension.
    :param ws_mult: List of multipliers (of #Walks) to be used for window size.
    :param evaluation: Column used to evaluate the hyperparameters ('auc', 'aps', ...)
    :return: The index to the parameter with the best result.
    """
    best_idx = -1
    best_val = 0.0

    if len([x for x in [nwal, wlen, dim, ws_mult] if len(x) != 1]) != 1:
        logger.error('One and only one parameter should be optimized at a time')
        return -1
    df = pd.DataFrame()
    fig, axs = plt.subplots()
    opt_start = time.time()
    for idx, (nw, wl, d, ws) in enumerate(itertools.product(nwal, wlen, dim, ws_mult)):
        auc_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            num_walks=nw,
            walk_length=wl,
            dimension=d,
            window_size=int(ceil(ws * nw)),
        )
        logger.debug(f'Intermediate result in {time.time() - opt_start}s')
        iter_string = f'nw{nw}wl{wl}d{d}ws{ws}'
        df[iter_string] = auc_df[evaluation]
        if auc_df[evaluation].mean() >= best_val:
            best_idx = idx
            best_val = auc_df[evaluation].mean()
    logger.info(f'Optimization {fig.number} total time: {time.time() - opt_start}s')
    df.boxplot(return_type='axes', ax=axs)
    fig.suptitle(f'Optimization #{fig.number}')
    fig.savefig(f'optimization_ad_{fig.number}.png')
    return best_idx


def main():

    dim = len(graph_paths)
    fig, axs = plt.subplots(ncols=dim, sharey='all', squeeze=False)
    fig.set_size_inches(15, 5)

    gene_list = parse_dge(
        dge_path=dge_path,
        entrez_id_header=entrez_id_name_ad,
        log2_fold_change_header=log_fold_change_name_ad,
        adj_p_header=adjusted_p_value_name_ad,
        entrez_delimiter=split_char,
        base_mean_header=base_mean_name_ad,
    )
    for ppi_graph_path in graph_paths:

        max_log2_fold_change, min_log2_fold_change = lfc_cutoff, lfc_cutoff * -1

        fig.suptitle(f'Hyperparameter optimization for embedding of {ppi_graph_path} using DGE {dge} ')
        network = generate_ppi_network(
            ppi_graph_path=os.path.join(ppi_base_path, ppi_graph_path),
            dge_list=gene_list,
            max_adj_p=max_padj,
            max_log2_fold_change=max_log2_fold_change,
            min_log2_fold_change=min_log2_fold_change,
            ppi_edge_min_confidence=ppi_edge_min_confidence,
            current_disease_ids_path='',
            disease_associations_path=phewas_path,
        )

        targets = parse_gene_list(targets_file, network.graph)

        # File with no weights
        write_gat2vec_input_files(
            network=network,
            targets=targets,
            home_dir=g2v_path
        )

        # Starting parameters
        n_walk = 10
        w_size = .125
        dimens = 128

        # Optimize walk length
        logger.info('Optimizing walk length')
        best_idx = optimize_parameters(
            network=network,
            nwal=[n_walk],
            wlen=g2v_opt_walk_len,
            dim=[dimens],
            ws_mult=[w_size],
            evaluation='auc'
        )
        w_len = g2v_opt_walk_len[best_idx]

        # Optimize number of walks
        logger.info('Optimizing number of walks')
        best_idx = optimize_parameters(
            network=network,
            nwal=g2v_opt_num_walks,
            wlen=[w_len],
            dim=[dimens],
            ws_mult=[w_size],
            evaluation='auc'
        )
        n_walk = g2v_opt_num_walks[best_idx]

        # Optimize window size
        logger.info('Optimizing window size')
        best_idx = optimize_parameters(
            network=network,
            nwal=[n_walk],
            wlen=[w_len],
            dim=[dimens],
            ws_mult=g2v_opt_win_size_mult,
            evaluation='auc'
        )
        w_size = g2v_opt_win_size_mult[best_idx]

        # Optimize dimension
        logger.info('Optimizing dimension')
        best_idx = optimize_parameters(
            network=network,
            nwal=[n_walk],
            wlen=[w_len],
            dim=g2v_opt_dimension,
            ws_mult=[w_size],
            evaluation='auc'
        )
        dimens = g2v_opt_dimension[best_idx]

        logger.info(f'Best Parameters for {ppi_graph_path}:')
        logger.info(f'walk length: {w_len}')
        logger.info(f'Number of walks: {n_walk}')
        logger.info(f'window size: {w_size}')
        logger.info(f'dimension: {dimens}')

        auc_df, _ = rank_targets(
            directory=g2v_path,
            network=network,
            num_walks=n_walk,
            walk_length=w_len,
            dimension=dimens,
            window_size=w_size,
        )

        logger.info('Result with optimized parameters.')
        logger.info(auc_df)

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
