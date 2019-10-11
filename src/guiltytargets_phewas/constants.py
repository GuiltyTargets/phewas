# -*- coding: utf-8 -*-

"""Constants for GuiltyTargets PheWAS."""

DATA_BASE_DIR = r'C:/users/mauricio/thesis/data'
AD_DGE_DATASETS = ['RosMap', 'Mayo_CBE', 'Mayo_CTX', 'BM10', 'BM22', 'BM36', 'BM44']
NON_AD_DGE_DATASETS = ['hc', 'lc', 'aml', 'ipf', 'ms']
DGE_DATASETS = {
    'ad': AD_DGE_DATASETS,
    'nonad': NON_AD_DGE_DATASETS
}
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
dge_params_ad = {
    'mean': base_mean_name_ad,
    'l2f': log_fold_change_name_ad,
    'adjp': adjusted_p_value_name_ad,
    'id': entrez_id_name_ad
}

# for differential expression (other files)
base_mean_name_dis = None
log_fold_change_name_dis = 'logFC'
adjusted_p_value_name_dis = 'adj.P.Val'
entrez_id_name_dis = 'Gene.ID'
dge_params_dis = {
    'mean': base_mean_name_dis,
    'l2f': log_fold_change_name_dis,
    'adjp': adjusted_p_value_name_dis,
    'id': entrez_id_name_dis
}

# for differential expression (any)
max_padj = 0.05
split_char = '///'
diff_type = 'all'

# STRING database
string_host = "localhost"
string_database = "stringitems"
string_user = "postgres"
string_password = "123456"
