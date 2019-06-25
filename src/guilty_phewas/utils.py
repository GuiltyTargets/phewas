# -*- coding: utf-8 -*-

"""Utility functions to run the algorithm"""

import os
from typing import Dict, List, Optional

import mygene
from networkx import DiGraph
from opentargets import OpenTargetsClient
import pandas as pd
from pybel.dsl import gene, protein, rna, BaseEntity

from .constants import disease_ids_efo, disease_abr, data_dir, ot_file


def get_base_dir(basedir, path):
    return os.path.join(basedir, path)


def get_local_bel_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return os.path.join(dir_, path + '.bel')


def get_struct_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return os.path.join(dir_, path + '_graph.adjlist')


def get_attr_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return os.path.join(dir_, path + '_na.adjlist')


def get_labels_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return os.path.join(dir_, 'labels_maped.txt')


def add_disease_attribute(graph: DiGraph, att_mappings: Dict):
    """Add the phenotypes to the Base Entities as attributes."""
    for node in graph.nodes:
        if ((isinstance(node, protein) or
             isinstance(node, rna) or
             isinstance(node, gene)) and
                node.name in att_mappings):
            graph.nodes[node]['phenotypes'] = [phtype for _, phtype in att_mappings[node.name]]


def write_adj_file_attribute(graph, filepath: str, att_mappings: Dict, pred_mapping: Optional[Dict] = None):
    """Write an adjacency file from the attribute graph."""
    if pred_mapping is None:
        pred_mapping = dict()
    with open(filepath, 'w') as f:
        for node in graph.nodes:
            line = f"{node}"
            if 'phenotypes' in graph.nodes[node]:  # "There are diseases in the node":
                line += f" {' '.join(str(att_mappings[phe]) for phe in graph.nodes[node]['phenotypes'])}"
            if 'label' in graph.nodes[node] and graph.nodes[node]['label'] in pred_mapping:
                line += f" {pred_mapping[graph.nodes[node]['label']]}"
            print(line, file=f)


# Copied from GuiltyTargets/reproduction
def download_for_disease(disease_id, outpath):
    ot = OpenTargetsClient()
    assoc = ot.get_associations_for_disease(
        disease_id,
        fields=['associationscore.datatypes', 'target.id']
    ).filter(
        datatype='known_drug'
    )
    ensembl_list = [a['target']['id'] for a in assoc]

    mg = mygene.MyGeneInfo()
    id_mappings = mg.getgenes(ensembl_list, fields="symbol")

    with open(outpath, 'w+') as outfile:
        for mapping in id_mappings:
            if 'symbol' in mapping.keys():
                outfile.write(mapping['symbol'])
                outfile.write('\n')


def parse_gene_list(path: str, graph: DiGraph, anno_type: str = "name") -> list:
    """Parse a list of genes and return them if they are in the network.

    :param str path: The path of input file.
    :param Graph graph: The graph with genes as nodes.
    :param str anno_type: The type of annotation with two options:name-Entrez ID, symbol-HGNC symbol.
    :return list: A list of genes, all of which are in the network.
    """
    # read the file
    genes = pd.read_csv(path, header=None)[0].tolist()
    genes = set(genes)

    # get those genes which are in the network
    symbols = [node['name'] for node in graph.nodes if isinstance(node, BaseEntity) and 'name' in node]

    return list(genes.intersection(symbols))


def write_labels_file(outfile: str, graph: DiGraph, target_list: List[str]):
    """Write the file with the nodes' ids and the labels (known target vs not).

    :param outfile: Path to the output file.
    :param graph: Graph (with the nodes already converted to integers).
    :param target_list: List with the known targets symbols.
    :return:
    """
    with open(outfile, 'w') as file:
        for node in graph.nodes:
            if (
                isinstance(graph.nodes[node]['old_label'], BaseEntity) and
                'name' in graph.nodes[node]['old_label'] and
                graph.nodes[node]['old_label']['name'] in target_list
            ):
                print(f'{node}\t1', file=file)
            else:
                print(f'{node}\t0', file=file)


if __name__ == '__main__':
    for (id, abr) in zip(disease_ids_efo, disease_abr):
        outpath = os.path.join(data_dir, abr, ot_file)
        download_for_disease(id, outpath)
