# -*- coding: utf-8 -*-

"""Runs the link prediction analysis to assess new disease-target associations."""

# Part5

from collections import defaultdict
from copy import deepcopy
import itertools as itt
import logging
import os
from time import time
from typing import List, Tuple

from GAT2VEC.parsers import get_embeddingDF
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from guiltytargets.gat2vec import Gat2Vec
from guiltytargets.constants import gat2vec_config
from guiltytargets.ppi_network_annotation import parse_dge
from guiltytargets_phewas.constants import *
from guiltytargets_phewas.utils import timed_main_run
from guiltytargets_phewas.target_repositioning import calculate_prob, generate_heterogeneous_network

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('link_prediction_cv.log')
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
g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'linkprediction3')
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


def extract_results(results_dict, lp_results, dataset, evaluation):
    for i, (auc, aps) in enumerate(lp_results):
        results_dict['tr'].append(i)
        results_dict['auc'].append(auc)
        results_dict['aps'].append(aps)
        results_dict['dge'].append(dataset)
        results_dict['eval'].append(evaluation)


def predict_links_cv(
        network,
        dataset,
        num_walks,
        walk_length,
        dimension,
        window_size
) -> pd.DataFrame:
    results = defaultdict(list)
    disease_abv = dataset_to_disease_abv(dataset)
    do_id = disease_identifiers[disease_abv]
    for i in range(10):
        for cv_labels in network.write_gat2vec_cv_split(
                home_dir=g2v_path,
                disease_id=do_id,
                filter_pleiotropic_targets=True
        ):
            g2v = Gat2Vec(
                input_dir=g2v_path,
                output_dir=g2v_path,
                label=False,
                tr=gat2vec_config.training_ratio
            )
            model = g2v.train_gat2vec(
                num_walks,
                walk_length,
                dimension,
                window_size,
                output=True,
            )
            model_df = get_embeddingDF(model)
            labels = pd.DataFrame(
                data=[
                    (x, y)
                    for x, y
                    in cv_labels.items()
                ],
                columns=['gene', 'label']
            )
            disease_idx = network.get_index_for_disease(do_id)
            labels['prob'] = pd.Series([
                calculate_prob(model_df[disease_idx], model_df[gene_idx])
                for gene_idx
                in labels.loc[:, 'gene']
            ])
            auc = roc_auc_score(labels['label'], labels['prob'])
            aps = average_precision_score(labels['label'], labels['prob'])
            results['TR'].append(i)
            results['dge'].append(dataset)
            results['eval'].append('5fold')
            results['auc'].append(auc)
            results['aps'].append(aps)

    df = pd.DataFrame(results)
    return df.groupby(axis=0, by="TR").mean()


def main():
    # natural order: disease <-> target <-> chem
    # disease - chem is what is desired
    # disease - target is what is desired
    # http://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz

    results = pd.DataFrame()
    h_network = generate_heterogeneous_network(
        ppi_path,
        dg_path,
        dch_path,
        chg_file
    )

    for use_dge, dataset in itt.product([False, True], ['Mayo_CTX']):
        disease_abv = dataset_to_disease_abv(dataset)
        dge_params = dge_params_ad if disease_abv == 'ad' else dge_params_dis
        do_id = disease_identifiers[disease_abv]
        logger.debug(f'Running for disease {disease_abv}, with the dataset {dataset}, '
                     f'using the id {disease_identifiers[disease_abv]}')
        try:
            gene_list = parse_dge(
                dge_path=dge_file(dataset),
                entrez_id_header=dge_params['id'],
                log2_fold_change_header=dge_params['l2f'],
                adj_p_header=dge_params['adjp'],
                entrez_delimiter=split_char,
                base_mean_header=dge_params['mean'],
            )

            h_network_ds = deepcopy(h_network)
            h_network_ds.set_up_network(genes=gene_list)

            num_walks = gat2vec_config.num_walks
            walk_length = gat2vec_config.walk_length
            dimension = gat2vec_config.dimension
            window_size = gat2vec_config.window_size

            metrics_df = predict_links_cv(
                h_network_ds,
                dataset=dataset,
                num_walks=num_walks,
                walk_length=walk_length,
                dimension=dimension,
                window_size=window_size
            )
            results = results.append(metrics_df, ignore_index=True)
        except ValueError:
            logger.error(f'Dataset {dataset} ({do_id}) not found in the graph.')

    results.to_csv(os.path.join(g2v_path, 'results_df.tsv'), sep='\t')

    print("done")


if __name__ == "__main__":
    timed_main_run(main, logger)
