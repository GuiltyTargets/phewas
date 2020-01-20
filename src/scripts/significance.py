# -*- coding: utf-8 -*-

"""Gets the results produced by other scripts for analysis."""

import logging
import os

import pandas as pd
from scipy.stats import f_oneway, ttest_ind
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.multitest import multipletests

from guiltytargets_phewas.constants import AD_DGE_DATASETS, DATA_BASE_DIR, NON_AD_DGE_DATASETS

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('anova.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

resuls_base_folder = os.path.join(DATA_BASE_DIR, 'gat2vec_files')

resuls_base_folder = r'c:/users/mauricio/thesis/final results'

results_subfolders = [
    'part0',
    'part1',
    'part2',
    'part3',
    'part4',
    'part5',
    'part6',
]

results_subfolders = [
    'part2',
    'part3',
]

"""lacerda@abi-srv-03:~/data/gat2vec_files/part2$ ls
labels_maped.txt   part2_graph.adjlist  results_df_ad.tsv
part2_gat2vec.emb  part2_na.adjlist     results_df_nonad.tsv"""


def get_dataframe(from_path: str) -> pd.DataFrame:
    if os.path.isfile(os.path.join(from_path, 'results_df.tsv')):
        return pd.read_csv(os.path.join(from_path, 'results_df.tsv'), sep='\t')
    else:
        df_ad = pd.read_csv(os.path.join(from_path, 'results_df_ad.tsv'), sep='\t')
        df_nonad = pd.read_csv(os.path.join(from_path, 'results_df_nonad.tsv'), sep='\t')
        return df_ad.append(df_nonad)


def two_way_anova(df: pd.DataFrame, file):
    """Calculates the statistical significance of the compared hyperparameter/ model being evaluated in each
    part (eval) and the differential gene expression/disease using a two-way ANOVA.
    https://pythonfordatascience.org/anova-2-way-n-way/

    :param df:
    :return:
    """
    for metric in ['auc', 'aps']:
        type1 = f'{metric} ~ C(eval) + C(dge) + C(eval)*C(dge)'
        model = ols(type1, data=df).fit()
        print(
            f"Overall model for {metric} F({model.df_model: .0f},{model.df_resid: .0f}) "
            f"= {model.fvalue: .3f}, p = {model.f_pvalue: .4f}",
            file=file
        )
        # Creates the ANOVA table
        res = sm.stats.anova_lm(model, typ=1)
        print(res, file=file)


def serial_ttest(df, eval_filter):
    ttests, pvals = list(), list()
    for dge in set(df['dge']):
        df_dge = df[df['dge'] == dge]
        t, p = ttest_ind(
            df_dge[df_dge['eval'] == eval_filter]['auc'],
            df_dge[df_dge['eval'] != eval_filter]['auc'],
            equal_var=False
        )
        ttests.append((dge, t))
        pvals.append(p)
    _, adjpvals, _, _ = multipletests(pvals, alpha=0.025, method='bonferroni')
    g1, g2 = [], []
    pval1, pval2 = 0., 0.
    for adjp, (d, t) in zip(adjpvals, ttests):
        if adjp < 0.025:
            if t > 0:
                g1.append((d, t))
                if pval1 < adjp:
                    pval1 = adjp
            else:
                g2.append((d, t))
                if pval2 < adjp:
                    pval2 = adjp
    return g1, pval1, g2, pval2


def main():
    """"""
    with open('anova.txt', 'w') as file:
        #############
        # NETWORK
        subfolder = 'part1'
        print(f'Processing {subfolder}:', file=file)
        df1 = get_dataframe(os.path.join(resuls_base_folder, subfolder))

        print(f'AD only', file=file)
        df = df1[df1['dge'].isin(AD_DGE_DATASETS)]
        two_way_anova(df, file)
        mc = MultiComparison(df['auc'], df['eval'])
        res1 = mc.tukeyhsd()
        print(res1.summary(), file=file)

        print(f'NONAD only', file=file)
        df = df1[df1['dge'].isin(NON_AD_DGE_DATASETS)]
        two_way_anova(df, file)
        g1, pval1, g2, pval2 = serial_ttest(df, 'st')
        print(f'g1 greater {g1}, p<{pval1}', file=file)
        print(f'g2 greater {g2}, p<{pval2}', file=file)

        #############
        # Classification
        subfolder = 'part2'
        print(f'Processing {subfolder}:', file=file)
        df = get_dataframe(os.path.join(resuls_base_folder, subfolder))
        two_way_anova(df, file)
        print('AUC-ROC', file=file)
        mc = MultiComparison(df['auc'], df['eval'])
        res1 = mc.tukeyhsd()
        print(res1.summary(), file=file)
        print('AUC-PR', file=file)
        mc = MultiComparison(df['aps'], df['eval'])
        res1 = mc.tukeyhsd()
        print(res1.summary(), file=file)

        #############
        # Weighting
        subfolder = 'part3'
        print(f'Processing {subfolder}:', file=file)
        df = get_dataframe(os.path.join(resuls_base_folder, subfolder))
        two_way_anova(df, file)
        g1, pval1, g2, pval2 = serial_ttest(df, 'weighted')
        print(f'g1 greater {g1}, p<{pval1}', file=file)
        print(f'g2 greater {g2}, p<{pval2}', file=file)

        #############
        # PheWAS
        subfolder = 'part6'
        print(f'Processing {subfolder} (PheWAS):', file=file)
        df = get_dataframe(os.path.join(resuls_base_folder, subfolder))
        two_way_anova(df, file)
        g1, pval1, g2, pval2 = serial_ttest(df, 'phewas')
        print(f'g1 greater {g1}, p<{pval1}', file=file)
        print(f'g2 greater {g2}, p<{pval2}', file=file)

        #############
        # link prediction
        subfolder = 'link_prediction'
        print(f'Processing {subfolder}:', file=file)
        df = get_dataframe(os.path.join(resuls_base_folder, subfolder))
        two_way_anova(df, file)
        g1, pval1, g2, pval2 = serial_ttest(df, True)
        print(f'g1 greater {g1}, p<{pval1}', file=file)
        print(f'g2 greater {g2}, p<{pval2}', file=file)

        #############
        # gat2vec
        subfolder = 'part4'
        print(f'Processing {subfolder}:', file=file)
        df2 = get_dataframe(os.path.join(resuls_base_folder, subfolder))
        for param in ['Dimension', 'Num Walks', 'Walk Length', 'Window Size']:
            print(f'{param} only', file=file)
            df = df2[df2['param'] == param]
            two_way_anova(df, file)
            print('AUC-ROC', file=file)
            mc = MultiComparison(df['auc'], df['eval'])
            res1 = mc.tukeyhsd()
            print(res1.summary(), file=file)
            print('AUC-PR', file=file)
            mc = MultiComparison(df['aps'], df['eval'])
            res1 = mc.tukeyhsd()
            print(res1.summary(), file=file)


if __name__ == '__main__':
    main()
