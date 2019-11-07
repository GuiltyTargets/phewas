# -*- coding: utf-8 -*-

"""Gets the results produced by other scripts for analysis."""

import os

import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp

from guiltytargets_phewas.constants import DATA_BASE_DIR

resuls_base_folder = os.path.join(DATA_BASE_DIR, 'gat2vec_files')

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
]


def i_suck(df):
    """https://pythonfordatascience.org/anova-2-way-n-way/"""
    type1 = 'auc ~ C(eval) + C(dge) + C(eval)*C(dge)'
    model = ols(type1, data=df).fit()
    print(
        f"Overall model F({model.df_model: .0f},{model.df_resid: .0f}) = {model.fvalue: .3f}, p = {model.f_pvalue: .4f}"
    )

    # Creates the ANOVA table
    res = sm.stats.anova_lm(model, typ=1)
    print(res)


def main():
    """"""
    for subfolder in results_subfolders:
        results_file = os.path.join(resuls_base_folder, subfolder, 'results_df.tsv')
        df = pd.read_csv(results_file, sep='\t')
        i_suck(df)

    # for each part:
    # 1) read the results_df.
    # 2) run the f-test.
    # 3) save the f-test result.


if __name__ == '__main__':
    main()
