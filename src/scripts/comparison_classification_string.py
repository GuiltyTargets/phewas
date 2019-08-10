# -*- coding: utf-8 -*-

"""Runs the ranking of drug targets using the STRING data as PPI network. It aims to compare the different
classification methods applied in this work."""

import logging
import os
import traceback
import warnings

import matplotlib.pyplot as plt
import pandas as pd

from guiltytargets_phewas.utils import download_for_disease, get_association_scores
from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_association_scores, parse_gene_list

# Suppress warnings
warnings.simplefilter('ignore')

# Log the run
logger = logging.getLogger('guiltytargets')
fh = logging.FileHandler('comp_class_string.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
logger.addHandler(ch)
logger.addHandler(fh)

# Paths
data_basedir = r'/home/bit/lacerda/data'
base_ppi_path = os.path.join(data_basedir, 'STRING')

targets_file_mask = os.path.join(data_basedir, r'OpenTargets/%s/ot_symbol.txt')
assoc_file_mask = os.path.join(data_basedir, r'OpenTargets/%s/ot_assoc_symbol.txt')
g2v_path_mask = os.path.join(data_basedir, r'gat2vec_files/ppi/%s')
phewas_path = r'phewas_catalog/phewas.txt'

dge_path_mask = os.path.join(data_basedir, r'DGE/%s/DifferentialExpression.tsv')

graph_paths = ['string_symbol.edgelist']

DISEASE_ABBREVIATIONS = [
    'aml',  # acute myeloid leukemia # bone marrow
    'hc',  # hepatocellular carcinoma # liver
    'ipf',  # idiopathic pulmonary fibrosis # lung
    'lc',  # cirrhosis of liver # liver
    'ms',  # multiple sclerosis # nervous system?
]
EFOS = {
    'aml': 'EFO_0000222',
    'hc': 'EFO_0000182',
    'ipf': 'EFO_0000768',
    'lc': 'EFO_0001422',
    'ms': 'EFO_0003885',
}
lfc_cutoff = 1.5  # no significance when changed
ppi_edge_min_confidence = 0.63  # confidence_cutoffs = [0.0, 0.63]

# for differential expression
max_padj = 0.05
base_mean_name = None
log_fold_change_name = 'logFC'
adjusted_p_value_name = 'adj.P.Val'
entrez_id_name = 'Gene.symbol'
split_char = '///'
diff_type = 'all'


def main():
    for disease in DISEASE_ABBREVIATIONS:
        targets_file = targets_file_mask % disease
        assoc_file = assoc_file_mask % disease
        g2v_path = g2v_path_mask % disease
        dge_path = dge_path_mask % disease
        print(f'downloading for {disease}')
        
        download_for_disease(EFOS[disease], targets_file, anno_type='symbol')
        get_association_scores(EFOS[disease], assoc_file, anno_type='symbol')
        print(f'downloaded for {disease}')

        gene_list = parse_dge(
            dge_path=dge_path,
            entrez_id_header=entrez_id_name,
            log2_fold_change_header=log_fold_change_name,
            adj_p_header=adjusted_p_value_name,
            entrez_delimiter=split_char,
            base_mean_header=base_mean_name,
        )

        dim = len(graph_paths)
        fig, axs = plt.subplots(ncols=dim, sharey='all')
        fig.set_size_inches(15, 5)
        fig.suptitle(f'DGE {disease}')

        df = pd.DataFrame()

        axs_ind = 0
        for ppi_graph_path in graph_paths:
            max_log2_fold_change, min_log2_fold_change = lfc_cutoff, lfc_cutoff * -1

            network = generate_ppi_network(
                ppi_graph_path=os.path.join(base_ppi_path, ppi_graph_path),
                dge_list=gene_list,
                max_adj_p=max_padj,
                max_log2_fold_change=max_log2_fold_change,
                min_log2_fold_change=min_log2_fold_change,
                ppi_edge_min_confidence=ppi_edge_min_confidence,
                current_disease_ids_path='',
                disease_associations_path=phewas_path,
            )

            targets = parse_gene_list(targets_file, network.graph)

            assoc_score = assoc_file and parse_association_scores(assoc_file)

            # File with no weights
            write_gat2vec_input_files(
                network=network,
                targets=targets,
                home_dir=g2v_path
            )

            auc_df, _ = rank_targets(
                directory=g2v_path,
                network=network,
            )

            df['rr'] = auc_df['auc']

            auc_df, _ = rank_targets(
                directory=g2v_path,
                network=network,
                evaluation='svm',
                class_weights='balanced'
            )

            df['bsvm'] = auc_df['auc']

            # File with weights
            write_gat2vec_input_files(
                network=network,
                targets=targets,
                home_dir=g2v_path,
                assoc_score=assoc_score
            )

            auc_df, _ = rank_targets(
                directory=g2v_path,
                network=network,
            )

            df['wrr'] = auc_df['auc']

            auc_df, _ = rank_targets(
                directory=g2v_path,
                network=network,
                evaluation='svm',
                class_weights='balanced'
            )

            df['wbsvm'] = auc_df['auc']

            df.boxplot(column=['rr', 'wrr', 'bsvm', 'wbsvm'], ax=axs)

            axs.set_title(f'PPI {ppi_graph_path}"')
            axs_ind += 1
        fig.savefig(f'comparison_humanbase-{disease}.png')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.error(type(e))
        logger.error(traceback.format_exc())
