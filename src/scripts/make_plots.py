# -*- coding: utf-8 -*-

"""Generate all plots from the saved dataframe with results."""


import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from guiltytargets_phewas.constants import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('plots.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.addHandler(fh)

results_base_folder = r'c:/users/mauricio/thesis/final results'

results_subfolders = [
    'part6',
    'part1',
    'part2',
    'part3',
    'part4',
    'part6',
]

assert os.path.isdir(results_base_folder)


change_labels = {
    'st': 'STRING',
    'bel': 'BEL',
    'hb': 'HumanBase',
    'ctex': 'C. Cortex',
    'hipp': 'Hippocampus',
    'hc': 'Hepatocellular\ncarcinoma',
    'lc': 'Liver\ncirrhosis',
    'ipf': 'Idiopathic\npulmonary\nfibrosis',
    'ms': 'Multiple\nsclerosis',
    'aml': 'Acute\nmyeloid\nleukemia',
    'BM10': 'MSBB BM\narea 10',
    'BM22': 'MSBB BM\narea 22',
    'BM36': 'MSBB BM\narea 36',
    'BM44': 'MSBB BM\narea 44',
    'RosMap': 'RosMap',
    'Mayo_CTX': 'Mayo\nTemporal\nCortex',
    'Mayo_CBE': 'Mayo\nCerebellum',
}


def dataset_to_disease_abv(dataset: str) -> str:
    return dataset if dataset in NON_AD_DGE_DATASETS else 'ad'


def get_dataframe(from_path: str) -> pd.DataFrame:
    if os.path.isfile(os.path.join(from_path, 'results_df.tsv')):
        df = pd.read_csv(os.path.join(from_path, 'results_df.tsv'), sep='\t')
    else:
        df_ad = pd.read_csv(os.path.join(from_path, 'results_df_ad.tsv'), sep='\t')
        df_nonad = pd.read_csv(os.path.join(from_path, 'results_df_nonad.tsv'), sep='\t')
        df = df_ad.append(df_nonad)
    df['dge_relabel'] = [change_labels[x] for x in df['dge']]
    df['sort1'] = [
        dataset_to_disease_abv(x)
        for x in df['dge']
    ]
    return df


def handle_plot(bp, xloc, yloc, legend_title, use_ncol=False):
    bp.set_xlabel("")
    ylabel = bp.get_ylabel()
    bp.set_ylabel("AUC-ROC" if ylabel == 'auc' else "AUC-PR")
    handles, labels = bp.get_legend_handles_labels()
    bp.legend(
        handles,
        labels,
        ncol=len(labels) if use_ncol else 1,
        loc=(xloc, yloc),
        title=legend_title
    )
    bp.set_xticklabels(
        bp.get_xticklabels(),
        rotation=30
    )


def main():
    """"""
    thiner = 0.5
    subfolder = 'part1'
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results['eval'] = [change_labels[x] for x in results['eval']]
    results['sort2'] = [
        0 if x == 'STRING' else 2 if x == 'BEL' else 1
        for x in results['eval']
    ]
    results = results.sort_values(by=['sort1', 'sort2', 'dge_relabel'])
    # Prepare results
    for metric in ['auc', 'aps']:
        for key, datasets in DGE_DATASETS.items():
            bp = sns.boxplot(
                data=results[results['dge'].isin(datasets)],
                x='dge_relabel',
                y=metric,
                hue='eval',
                width=thiner if key == 'nonad' else 1
            )
            handle_plot(bp, 0.1, 0.45, 'Network')
            fig = bp.get_figure()
            fig.tight_layout()
            fig.savefig(os.path.join(results_base_folder, f'compare1_{metric}_{key}.png'))  # AUC non-AD and AD
            plt.close()
    subfolder = 'part3'  # Weighting
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results['eval'] = [
        'Yes' if x == 'weighted' else 'No'
        for x in results['eval']
    ]
    results = results.sort_values(by=['sort1', 'dge_relabel', 'eval'])
    # Prepare results
    for metric in ['auc', 'aps']:
        bp = sns.boxplot(
            data=results,
            x='dge_relabel',
            y=metric,
            hue='eval',
            width=thiner
        )
        handle_plot(bp, 0.05, 0.3, 'Weighting')
        fig = bp.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(left=0.07)
        fig.set_size_inches(12.8, 4.8)
        fig.savefig(os.path.join(results_base_folder, f'compare4_{metric}.png'))  # AUC non-AD and AD
        plt.close()
    subfolder = 'part4'  # g2v
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results = results.sort_values(by=['sort1', 'dge_relabel', 'eval'])
    # Prepare results
    for metric in ['auc', 'aps']:
        for analyzed_param in set(results['param']):
            colors = sns.color_palette()
            # swap blue to the default g2v param
            # dim = 128, nw = 10, wl = 80, ws = 5
            if analyzed_param in ['Dimension', 'Walk Length']:
                (colors[0], colors[2]) = (colors[2], colors[0])
            elif analyzed_param in ['Num Walks', 'Window Size']:
                (colors[0], colors[1]) = (colors[1], colors[0])
            bp = sns.boxplot(
                data=results[results['param'] == analyzed_param],
                x='dge_relabel',
                y=metric,
                hue='eval',
                palette=colors
            )
            handle_plot(bp, 0.01, 0.05, analyzed_param, use_ncol=True)
            fig = bp.get_figure()
            fig.tight_layout()
            fig.subplots_adjust(left=0.07)
            fig.set_size_inches(12.8, 4.8)
            param = ''.join(analyzed_param.split(' '))
            fig.savefig(os.path.join(results_base_folder, f'compare5_{metric}_{param}.png'))  # AUC non-AD and AD
            plt.close()
    subfolder = 'part6'
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results['eval'] = [
        'Yes' if x == 'phewas' else 'No'
        for x in results['eval']
    ]
    results['dge_relabel'] = [change_labels[x] for x in results['dge']]
    results = results.sort_values(by=['sort1', 'dge_relabel', 'eval'])
    # Prepare results
    for metric in ['auc', 'aps']:
        bp = sns.boxplot(
            data=results,
            x='dge_relabel',
            y=metric,
            hue='eval',
            width=thiner
        )
        handle_plot(bp, 0.1, 0.3, 'Using PheWAS')
        fig = bp.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(left=0.07)
        fig.set_size_inches(12.8, 4.8)
        fig.savefig(os.path.join(results_base_folder, f'compare2_{metric}.png'))  # AUC non-AD and AD
        plt.close()

    subfolder = 'part2'
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results['sort2'] = [
        0 if x == 'cv' else 1 if x == 'nested_cv' else 2
        for x in results['eval']
    ]
    results['eval'] = [
        'Logistic regression' if x == 'cv' else 'Nested logistic regression' if x == 'nested_cv' else 'Biased SVM'
        for x in results['eval']
    ]
    results['dge_relabel'] = [change_labels[x] for x in results['dge']]
    results = results.sort_values(by=['sort1', 'dge_relabel', 'sort2'])
    # Prepare results
    for metric in ['auc', 'aps']:
        bp = sns.boxplot(
            data=results,
            x='dge_relabel',
            y=metric,
            hue='eval',
            width=thiner + 0.1
        )
        handle_plot(bp, 0.01, 0.01, 'Classification')
        fig = bp.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(left=0.07)
        fig.set_size_inches(12.8, 4.8)
        fig.savefig(os.path.join(results_base_folder, f'compare3_{metric}.png'))  # AUC non-AD and AD
        plt.close()

    subfolder = 'link_prediction'  # Link prediction
    results = get_dataframe(os.path.join(results_base_folder, subfolder))
    results['eval'] = [
        'Yes' if x else 'No'
        for x in results['eval']
    ]
    results = results.sort_values(by=['sort1', 'dge_relabel', 'eval'])
    # Prepare results
    for metric in ['auc', 'aps']:
        bp = sns.boxplot(
            data=results,
            x='dge_relabel',
            y=metric,
            hue='eval',
            width=thiner
        )
        handle_plot(bp, 0.05, 0.3, 'Use DGE')
        fig = bp.get_figure()
        fig.tight_layout()
        fig.subplots_adjust(left=0.07)
        fig.set_size_inches(12.8, 4.8)
        fig.savefig(os.path.join(results_base_folder, f'compare7_{metric}.png'))  # AUC non-AD and AD
        plt.close()


if __name__ == '__main__':
    main()
