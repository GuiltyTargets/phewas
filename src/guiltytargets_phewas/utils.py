# -*- coding: utf-8 -*-

"""Utility functions to run the algorithm."""

import os
from typing import Dict, List, Optional

import bio2bel_phewascatalog
import mygene
import networkx as nx
import pandas as pd
from opentargets import OpenTargetsClient
from pybel.dsl import BaseEntity, gene, protein, rna

from guiltytargets.ppi_network_annotation.model.gene import Gene
from .network_phewas import NetworkNx


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


def add_disease_attribute(graph: nx.Graph, att_mappings: Dict):
    """Add the phenotypes to the Base Entities as attributes."""
    for node in graph:
        if isinstance(node, (protein, rna, gene) and node.name in att_mappings):
            graph.nodes[node]['phenotypes'] = [phtype for _, phtype in att_mappings[node.name]]


def add_dge_attribute(graph: nx.Graph, gene_dict: Dict):
    """ """
    max_adjp = 0.05
    min_l2fc, max_l2fc = -0.5, 0.5
    dge = {g.entrez_id or g.symbol: g.log2_fold_change for g in gene_dict if g.padj < max_adjp}

    for node in graph:
        if isinstance(node, (protein, rna, gene)) and node['name'] in dge.keys():
            if dge[node['name']] > max_l2fc:
                graph.nodes[node]['deg'] = True
                graph.nodes[node]['upreg'] = True
            elif dge[node['name']] < min_l2fc:
                graph.nodes[node]['deg'] = True
                graph.nodes[node]['downreg'] = True


def get_significantly_differentiated(gene_list: List[Gene], max_adjp: float):
    """Returns a dictionary only with significantly differentially expressed genes from the gene list."""
    max_adjp = max_adjp or 0.05

    dge = {
        g.entrez_id or g.symbol: g.log2_fold_change
        for g in gene_list
        if g.padj < max_adjp
    }

    return {k: v for k, v in dge.items() if k}


def write_adj_file_attribute(graph, filepath: str, att_mappings: Dict, pred_mapping: Optional[Dict] = None):
    """Write an adjacency file from the attribute graph."""
    if pred_mapping is None:
        pred_mapping = dict()
    with open(filepath, 'w') as file:
        for node in graph.nodes:
            line = f"{node}"
            if 'phenotypes' in graph.nodes[node]:  # There are diseases in the node (gene)
                line += f" {' '.join(str(att_mappings[phe]) for phe in graph.nodes[node]['phenotypes'])}"
            if 'label' in graph.nodes[node] and graph.nodes[node]['label'] in pred_mapping:  # predicate node
                line += f" {pred_mapping[graph.nodes[node]['label']]}"
            print(line, file=file)


# Copied from GuiltyTargets/reproduction
def download_for_disease(disease_id, outpath):
    ot = OpenTargetsClient()
    assoc = ot.get_associations_for_disease(
        disease_id,
        fields=['association_scoredatatypes', 'target.id']
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


def get_association_scores(disease_id, outpath):
    """Obtain the association scores from the specified disease that are
    stored in the OpenTargets database.

    :param disease_id: The EFO code to the disease.
    :param outpath: The path to the file to be created.
    :return:b
    """
    ot = OpenTargetsClient()
    assoc = ot.get_associations_for_disease(
        disease_id,
        fields=['association_scoreoverall', 'target.id']
    )
    assoc_simple = [
        {
            'id': a['target']['id'],
            'score': a['association_score']['overall']
        }
        for a in assoc
    ]
    ensembl_list = [a['id'] for a in assoc_simple]

    # Obtain the symbols for the genes associated to disease_id
    mg = mygene.MyGeneInfo()
    gene_query = mg.getgenes(ensembl_list, fields="symbol,query")
    # TODO change here for entrez id
    id_mappings = {
        g['query']: g['symbol']
        for g in gene_query
        if 'symbol' in g
    }
    # Get the symbols and the scores
    ensembl_list = [
        (id_mappings[a['id']], a['score'])
        for a in assoc_simple
        if a['id'] in id_mappings
    ]

    with open(outpath, 'w+') as outfile:
        for symbol, score in ensembl_list:
            print(f'{symbol}\t{score}', file=outfile)


def parse_gene_list(path: str, graph: nx.Graph, anno_type: str = "name") -> list:
    """Parse a list of genes and return them if they are in the network.

    :param str path: The path of input file.
    :param Graph graph: The graph with genes as nodes.
    :param str anno_type: The type of annotation with two options:name-Entrez ID, symbol-HGNC symbol.
    :return list: A list of genes, all of which are in the network.
    """
    # read the file
    df = pd.read_csv(
        path,
        names=['gene'],
        dtype={'gene': str}
    )
    genes = set(df['gene'])

    # get those genes which are in the network
    symbols = [
        node['name']
        for node in graph.nodes
        if isinstance(node, BaseEntity) and 'name' in node
    ]

    return list(genes.intersection(symbols))


def write_labels_file(outfile: str, graph: nx.Graph, target_list: List[str]):
    """Write the file with the nodes' ids and the labels (known target vs not).

    :param outfile: Path to the output file.
    :param graph: Graph (with the nodes already converted to integers).
    :param target_list: List with the known targets symbols.
    :return:
    """
    with open(outfile, 'w') as file:
        for node in graph.nodes:
            if (
                isinstance(graph.nodes[node]['bel_node'], BaseEntity) and
                'name' in graph.nodes[node]['bel_node'] and
                graph.nodes[node]['bel_node']['name'] in target_list
            ):
                print(f'{node}\t1', file=file)
            else:
                print(f'{node}\t0', file=file)


# parsers
def generate_phewas_file(out_path: str):
    """Create a file with gene-disease association from PheWAS data.

    :param out_path: Path where the file will be created.
    :return:
    """
    phewas_manager = bio2bel_phewascatalog.Manager()
    pw_dict = phewas_manager.to_dict()
    with open(out_path, 'w') as file:
        for target, diseaselist in pw_dict.items():
            for disease in diseaselist:
                print(f'{target} {disease[1]}', file=file)


# TODO delete this?
def refied_graph_to_ppi_adj_file(graph: nx.Graph, out_path: str):
    """Converts a Reified BEL Graph into a iGraph network, while keeping the connections and the gene information.

    :param graph: A reified BEL graph.
    :return: A iGraph version of the reified graph
    """
    converted = nx.Graph()
    for edge in graph.edges(key=True):
        (u, v, key) = edge
        print(u)
        print(v)
        print(key)
        print(edge)
    nx.write_adjlist(graph, out_path)


def get_phenotype_list(rbg) -> List:
    """ Get the list of phenotypes per gene in the graph.

    :param start_ind:
    :return:
    """
    pw_dict = bio2bel_phewascatalog.Manager().to_dict()
    start_ind = len(rbg.nodes) + 3
    phenotypes_list = set([phe for _list in pw_dict.values()
                           for odds, phe in _list])
    enum_phenotypes = enumerate(phenotypes_list, start=start_ind)
    phenot_hash = {phenot: num for num, phenot in enum_phenotypes}

    annotated_graph = rbg.copy()
    add_disease_attribute(annotated_graph, pw_dict)


class AttributeFileGenerator:

    def __init__(self, network: NetworkNx):
        """Initializes the object.

        :param graph: The graph with nodes already converted to integers.
        """
        self.graph = network.graph

    def get_attribute_mappings(self) -> Dict:
        mappings = dict()

        # Pega o numero de vertices do graph
        start_ind = len(self.graph)

        # att_mappings +1 +2 +3 vao pra DGE

        # as doencas comecam no +3

        # retorna dicionario q vai pro arquivo praticamente direto
        return mappings

    # TODO usar esse inves do write_adj_file_attribute
    def write_attribute_adj_list(self, outpath: str):
        """Write the bipartite attribute graph to a file.

        :param str outpath: Path to the output file.
        """
        att_mappings = self.get_attribute_mappings()

        with open(outpath, mode="w") as file:
            for k, v in att_mappings.items():
                print("{} {}".format(k, " ".join(str(e) for e in v)), file=file)

    def get_phenotype_list(self, rbg) -> List:
        """ Get the list of phenotypes per gene in the graph.

        :param start_ind:
        :return:
        """
        pw_dict = bio2bel_phewascatalog.Manager().to_dict()
        start_ind = len(rbg.nodes) + 3
        phenotypes_list = set([phe for _list in pw_dict.values()
                               for odds, phe in _list])
        enum_phenotypes = enumerate(phenotypes_list, start=start_ind)
        phenot_hash = {phenot: num for num, phenot in enum_phenotypes}

        annotated_graph = rbg.copy()
        add_disease_attribute(annotated_graph, pw_dict)

    @staticmethod
    def _add_attribute_values(value, att_mappings, indices):
        """Add an attribute value to the given vertices.

        :param int value: Attribute value.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :param list indices: Indices of the vertices.
        """
        for i in indices:
            att_mappings[i].append(value)
