# -*- coding: utf-8 -*-

"""Runs the link prediction analysis to assess new disease-target associations."""

# Part5

from collections import defaultdict
from copy import deepcopy
import itertools as itt
import logging
import multiprocessing as mp
import os
from time import time
from typing import List, Tuple

import pandas as pd

from guiltytargets.constants import gat2vec_config
from guiltytargets.ppi_network_annotation import parse_dge
from guiltytargets_phewas.constants import *
from guiltytargets_phewas.utils import timed_main_run
from guiltytargets_phewas.target_repositioning import generate_heterogeneous_network, predict_links

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('link_prediction2.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)
logger.addHandler(fh)

assert os.path.isdir(DATA_BASE_DIR), "Update your data_basedir folder for this environment."

# Paths
snap_path = os.path.join(DATA_BASE_DIR, 'SNAP')
chg_file = os.path.join(snap_path, 'ChG-Miner_miner-chem-gene.tsv.gz')
dch_path = os.path.join(snap_path, 'DCh-Miner_miner-disease-chemical.tsv.gz')
dg_path = os.path.join(snap_path, 'DG-AssocMiner_miner-disease-gene.tsv.gz')
ppi_path = os.path.join(DATA_BASE_DIR, 'STRING', 'string_entrez.edgelist')

targets_file = os.path.join(DATA_BASE_DIR, 'OpenTargets', 'ad', 'ot_symbol.txt')
g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'linkprediction2')
phewas_path = os.path.join(DATA_BASE_DIR, 'phewas_catalog', 'phewas_symbol.txt')

dge_base_path = os.path.join(DATA_BASE_DIR, 'DGE')


def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def dge_file(dge_code: str) -> str:
    file = 'DifferentialExpression' + ('.csv' if dataset_to_disease_abv(dge_code) == 'ad' else '.tsv')
    return os.path.join(dge_base_path, dge_code, file)


disease_identifiers = {
    'ad': 'DOID:10652',
    'lc': 'DOID:5082',
    'ipf': 'DOID:0050156',
    'ms': 'DOID:2377',
    'aml': 'DOID:9119',
    'hc': 'MESH:D006528',  # DOID:0070328, DOID:684 or DOID:5005
}


def mp_predict_links(
        num_walks: int,
        walk_length: int,
        dimension: int,
        window_size: int
) -> List[Tuple[float, float]]:
    pool = mp.Pool(mp.cpu_count())
    results_iter = [
        pool.apply(
            predict_links,
            args=(
                g2v_path,
                num_walks,
                walk_length,
                dimension,
                window_size
            )
        )
        for _
        in range(10)
    ]
    pool.close()
    pool.join()
    return results_iter


def extract_results(results_dict, lp_results, dataset, param, evaluation):
    for i, (auc, aps) in enumerate(lp_results):
        results_dict['tr'].append(i)
        results_dict['auc'].append(auc)
        results_dict['aps'].append(aps)
        results_dict['dge'].append(dataset)
        results_dict['eval'].append(evaluation)
        results_dict['param'].append(param)


def main():
    # natural order: disease <-> target <-> chem
    # disease - chem is what is desired
    # disease - target is what is desired
    # http://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz

    results_dict = defaultdict(list)
    h_network1 = generate_heterogeneous_network(
        ppi_path,
        dg_path,
        dch_path,
        chg_file
    )

    for use_dge, dataset in itt.product([True], AD_DGE_DATASETS + NON_AD_DGE_DATASETS):
        disease_abv = dataset_to_disease_abv(dataset)
        do_id = disease_identifiers[disease_abv]
        dge_params = dge_params_ad if disease_abv == 'ad' else dge_params_dis
        logger.debug(f'Running for disease {disease_abv}, with the dataset {dataset}, using the id {do_id}')
        try:
            gene_list = parse_dge(
                dge_path=dge_file(dataset),
                entrez_id_header=dge_params['id'],
                log2_fold_change_header=dge_params['l2f'],
                adj_p_header=dge_params['adjp'],
                entrez_delimiter=split_char,
                base_mean_header=dge_params['mean'],
            )

            h_network = deepcopy(h_network1)
            h_network.set_up_network(genes=gene_list)
            h_network.write_gat2vec_input_files(
                home_dir=g2v_path,
                disease_id=do_id,
                filter_pleiotropic_targets=True
            )

            # num_walks = gat2vec_config.num_walks
            walk_length = gat2vec_config.walk_length
            dimension = gat2vec_config.dimension
            window_size = gat2vec_config.window_size

            param = 'nw'
            for num_walks in [6, 10, 20, 40, 80]:
                start = time()
                lp_results = mp_predict_links(num_walks, walk_length, dimension, window_size)

                extract_results(results_dict, lp_results, dataset, param, num_walks)
                logger.info(f'Runtime for num_walks = {num_walks}: {time() - start}s')
            # best result from num_walks
            num_walks = gat2vec_config.num_walks

            param = 'wl'
            for walk_length in [20, 40, 80, 120, 160]:
                start = time()
                lp_results = mp_predict_links(num_walks, walk_length, dimension, window_size)

                extract_results(results_dict, lp_results, dataset, param, walk_length)
                logger.info(f'Runtime for walk_length = {walk_length}: {time() - start}s')
            # best result from num_walks
            walk_length = gat2vec_config.walk_length

            param = 'ws'
            for window_size in [3, 5, 7, 10, 20, 40]:
                start = time()
                lp_results = mp_predict_links(num_walks, walk_length, dimension, window_size)

                extract_results(results_dict, lp_results, dataset, param, window_size)
                logger.info(f'Runtime for window_size = {window_size}: {time() - start}s')
            # best result from num_walks
            window_size = gat2vec_config.window_size

            param = 'd'
            for dimension in [32, 64, 128, 256]:
                start = time()
                lp_results = mp_predict_links(num_walks, walk_length, dimension, window_size)

                extract_results(results_dict, lp_results, dataset, param, dimension)
                logger.info(f'Runtime for dimension = {dimension}: {time() - start}s')
        except ValueError:
            logger.error(f'Dataset {dataset} ({do_id}) not found in the graph.')

        results = pd.DataFrame(results_dict)
        # backup intermediate results in each iteration.
        results.to_csv(os.path.join(g2v_path, 'results_df_nonad.tsv'), sep='\t')

    print("done")


if __name__ == "__main__":
    timed_main_run(main, logger)
