# -*- coding: utf-8 -*-

"""Runs the ranking of drug targets using the HumanBase data as PPI network."""

import logging
import os
import time
import traceback
import warnings

import matplotlib.pyplot as plt
import pandas as pd

from guiltytargets.pipeline import rank_targets, write_gat2vec_input_files
from guiltytargets.ppi_network_annotation import generate_ppi_network, parse_dge
from guiltytargets.ppi_network_annotation.parsers import parse_association_scores, parse_gene_list

# Suppress warnings
warnings.simplefilter('ignore')

# Log the run
logger = logging.getLogger()
fh = logging.FileHandler('classify_diseases_hb.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

# Paths
data_base_dir = r'/home/bit/lacerda/data'
assert os.path.isdir(data_base_dir), "Update your data_basedir folder for this environment."
assert os.path.isdir(data_base_dir)

targets_file = os.path.join(data_base_dir, r'OpenTargets/ad/ot_entrez.txt')
assoc_file = os.path.join(data_base_dir, r'OpenTargets/ad/ot_assoc_entrez.txt')
g2v_path = os.path.join(data_base_dir, r'gat2vec_files/ppi/ad')
phewas_path = None  # Phewas file need to be converted to Entrez

ppi_base_path = os.path.join(data_base_dir, r'HumanBase')

dge_base_path = os.path.join(data_base_dir, r'DGE/AMP-AD/DifferentialExpression_%s.csv')

graph_paths = ['hippocampus_top', 'cerebral_cortex_top']

lfc_cutoff = 1.5  # no significance when changed
ppi_edge_min_confidence = 0.0

# for differential expression
max_padj = 0.05
base_mean_name = 'baseMean' or None  # it was None
log_fold_change_name = 'log2FoldChange'
adjusted_p_value_name = 'padj'
entrez_id_name = 'entrez'
split_char = '///'
diff_type = 'all'


def main():
    for dge in ['BM10', 'BM22', 'BM36', 'BM44']:

        dge_path = dge_base_path % dge

        gene_list = parse_dge(
            dge_path=dge_path,
            entrez_id_header=entrez_id_name,
            log2_fold_change_header=log_fold_change_name,
            adj_p_header=adjusted_p_value_name,
            entrez_delimiter=split_char,
            base_mean_header=base_mean_name,
            csv_separator=';'
        )

        dim = len(graph_paths)
        fig, axs = plt.subplots(ncols=dim, sharey='all', squeeze=False)
        fig.set_size_inches(10, 5)
        fig.suptitle(f'DGE {dge}')

        df = pd.DataFrame()

        axs_ind = 0
        for ppi_graph_path in graph_paths:
            max_log2_fold_change, min_log2_fold_change = lfc_cutoff, lfc_cutoff * -1

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

            df.boxplot(column=['rr', 'wrr', 'bsvm', 'wbsvm'], ax=axs[0][axs_ind])

            axs[0][axs_ind].set_title(f'PPI {ppi_graph_path}"')
            axs_ind += 1
        fig.savefig(f'comparison_humanbase({dge}).png')


if __name__ == '__main__':
    start_time = time.time()
    try:
        main()
    except Exception as e:
        logger.error(type(e))
        logger.error(traceback.format_exc())
    finally:
        logger.info(f"Total time: {time.time() - start_time}")
