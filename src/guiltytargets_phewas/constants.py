# -*- coding: utf-8 -*-

"""Constants for GuiltyTargets PheWAS."""

DATA_BASE_DIR = r'/home/bit/lacerda/data'
disease_abr = ['ad']
disease_ids_efo = ['EFO_0000249']
ot_file = 'ot_entrez.txt'

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
