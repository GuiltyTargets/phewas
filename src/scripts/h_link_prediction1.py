# -*- coding: utf-8 -*-

"""Runs the link prediction analysis to assess new disease-target associations."""

# Part8a

from collections import defaultdict
from copy import deepcopy
import itertools as itt
import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from guiltytargets.constants import gat2vec_config
from guiltytargets.ppi_network_annotation import parse_dge
from guiltytargets_phewas.constants import *
from guiltytargets_phewas.utils import timed_main_run
from guiltytargets_phewas.target_repositioning import generate_heterogeneous_network, predict_links

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('h_link_prediction.log')
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
disgenet_path = os.path.join(DATA_BASE_DIR, 'disgenet', 'curated_gene_disease_associations.tsv')
chg_file = os.path.join(snap_path, 'ChG-Miner_miner-chem-gene.tsv.gz')
dch_path = os.path.join(snap_path, 'DCh-Miner_miner-disease-chemical.tsv.gz')
dg_path = os.path.join(snap_path, 'DG-AssocMiner_miner-disease-gene.tsv.gz')
ppi_path = os.path.join(DATA_BASE_DIR, 'STRING', 'string_entrez.edgelist')

targets_file = os.path.join(DATA_BASE_DIR, 'OpenTargets', 'ad', 'ot_symbol.txt')
g2v_path = os.path.join(DATA_BASE_DIR, 'gat2vec_files', 'linkprediction')
phewas_path = os.path.join(DATA_BASE_DIR, 'phewas_catalog', 'phewas_symbol.txt')

dge_base_path = os.path.join(DATA_BASE_DIR, 'DGE')


def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def dge_file(dge_code: str) -> str:
    file = 'DifferentialExpression' + ('.csv' if dataset_to_disease_abv(dge_code) == 'ad' else '.tsv')
    return os.path.join(dge_base_path, dge_code, file)


disease_identifiers = {
    'ad': 'UMLS_CUI:C0750900',
    'lc': 'UMLS_CUI:C0023890',
    'ipf': 'UMLS_CUI:C1800706',
    'ms': 'UMLS_CUI:C0751324',
    'aml': 'UMLS_CUI:C0023467',
    'hc': 'UMLS_CUI:C2239176',  # DOID:0070328, DOID:684 or DOID:5005 or MESH:D006528
}


def main():

    results_dict = defaultdict(list)
    h_network1 = generate_heterogeneous_network(
        ppi_path,
        dg_path,
        dch_path,
        chg_file
    )

    for use_dge, dataset in itt.product([False, True], AD_DGE_DATASETS[2] + NON_AD_DGE_DATASETS):
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
                csv_separator=';'
            )

            h_network = deepcopy(h_network1)
            h_network.set_up_network(genes=gene_list)
            h_network.write_gat2vec_input_files(
                home_dir=g2v_path,
                disease_id=do_id,
                use_dge_data=use_dge,
                filter_pleiotropic_targets=True
            )
            for i in range(10):
                auc, aps = predict_links(
                    g2v_path,
                    gat2vec_config.num_walks,
                    gat2vec_config.walk_length,
                    gat2vec_config.dimension,
                    gat2vec_config.window_size,
                )
                logger.debug(f'tr: {i}\t{round(auc, 3)}\t{round(aps, 3)}')
                results_dict['tr'].append(i)
                results_dict['auc'].append(auc)
                results_dict['aps'].append(aps)
                results_dict['dge'].append(dataset)
                results_dict['eval'].append(use_dge)
        except ValueError:
            logger.error(f'Dataset {dataset} ({do_id}) not found in the graph.')
    results = pd.DataFrame(results_dict)
    results.to_csv(os.path.join(g2v_path, 'results_df.tsv'), sep='\t')

    for key, datasets in DGE_DATASETS.items():
        for metric in ['auc', 'aps']:
            fig = sns.boxplot(
                data=results[results['dge'].isin(datasets)],
                x='dge',
                y='auc',
                hue='eval'
            ).get_figure()
            fig.suptitle('Link Prediction.')
            fig.savefig(f'comp7_{key}_{metric}.png')
            plt.close()

        print("done")


if __name__ == "__main__":
    timed_main_run(main, logger)
