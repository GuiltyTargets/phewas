# -*- coding: utf-8 -*-

"""Utility functions to help link prediction."""


from math import exp
from typing import Tuple

from GAT2VEC.parsers import get_embeddingDF
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from guiltytargets.constants import gat2vec_config
from guiltytargets.gat2vec import Gat2Vec, gat2vec_paths
from guiltytargets.ppi_network_annotation.parsers import parse_ppi_graph
from guiltytargets_phewas.parsers import parse_disease_gene_graph, parse_disease_drug_graph, parse_gene_drug_graph
from guiltytargets_phewas.target_repositioning.heterogeneous_network import HeterogeneousNetwork


def calculate_prob(v1: pd.Series, v2: pd.Series) -> float:
    """Gets the p(u, v), given by the formula below."""
    return 1. / (1. + exp(- v1.dot(v2)))


def predict_links(
        input_path: str,
        num_walks: int,
        walk_length: int,
        dimension: int,
        window_size: int,
) -> Tuple[float, float]:
    g2v = Gat2Vec(input_path, input_path, label=False, tr=gat2vec_config.training_ratio)
    model = g2v.train_gat2vec(
        num_walks,
        walk_length,
        dimension,
        window_size,
        output=True,
    )
    df = get_embeddingDF(model)
    labels = pd.read_csv(
        gat2vec_paths.get_labels_path(input_path),
        sep='\t',
        header=None,
        names=['disease', 'gene', 'label']
    )
    disease_idx = labels.iloc[0, 0]
    labels['prob'] = pd.Series([
        calculate_prob(df[disease_idx], df[gene_idx])
        for gene_idx
        in labels.loc[:, 'gene']
    ])
    auc = roc_auc_score(labels['label'], labels['prob'])
    aps = average_precision_score(labels['label'], labels['prob'])
    return auc, aps


def generate_heterogeneous_network(
        ppi_graph_path,
        disease_gene_path: str,
        disease_drug_path: str,
        gene_drug_path: str,
) -> HeterogeneousNetwork:
    """"""
    disease_gene_interactions = parse_disease_gene_graph(disease_gene_path)
    disease_drug_interactions = parse_disease_drug_graph(disease_drug_path)
    gene_drug_interactions = parse_gene_drug_graph(gene_drug_path)

    protein_interactions = parse_ppi_graph(ppi_graph_path)
    protein_interactions = protein_interactions.simplify()

    network = HeterogeneousNetwork(
        protein_interactions,
        disease_gene_interactions,
        disease_drug_interactions,
        gene_drug_interactions
    )
    return network