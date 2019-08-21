# -*- coding: utf-8 -*-
import os

# Paths
DATA_BASE_DIR = r'/home/bit/lacerda/data'

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

# Values for hyper-parameter optimization
g2v_opt_num_walks = [2, 5, 10, 20, 30, 80]
g2v_opt_walk_len = [4, 20, 40, 80, 120, 160]
g2v_opt_dimension = [2, 16, 64, 128, 512]
g2v_opt_win_size = [2, 5, 10, 20, 30]

