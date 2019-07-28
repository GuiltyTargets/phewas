# -*- coding: utf-8 -*-

"""Utility functions to run the algorithm."""

from collections import defaultdict
import logging
from os import listdir
from os.path import isfile, join
from typing import Dict, List, Optional

import bio2bel_phewascatalog
import mygene
import networkx as nx
import pandas as pd
from opentargets import OpenTargetsClient
from pybel import BELGraph, from_path, union, collapse_to_genes
from pybel.dsl import BaseEntity, CentralDogma
from pybel_tools.assembler.reified_graph.assembler import reify_bel_graph

from guiltytargets.gat2vec import gat2vec_paths
from guiltytargets.ppi_network_annotation.model.gene import Gene
from guiltytargets.ppi_network_annotation.parsers import parse_disease_ids, parse_disease_associations
from .network_phewas import NetworkNx

logger = logging.getLogger('main')


def get_base_dir(basedir, path):
    return join(basedir, path)


def get_local_bel_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return join(dir_, path + '.bel')


def get_struct_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return join(dir_, path + '_graph.adjlist')


def get_attr_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return join(dir_, path + '_na.adjlist')


def get_labels_file(basedir, path):
    dir_ = get_base_dir(basedir, path)
    return join(dir_, 'labels_maped.txt')


def generate_bel_network(
        bel_graph_path: str,
        dge_list: List[Gene],
        max_adj_p: float,
        max_log2_fold_change: float,
        min_log2_fold_change: float,
        current_disease_ids_path: Optional[str] = None,
        disease_associations_path: Optional[str] = None,
) -> NetworkNx:
    """Generate the protein-protein interaction network.

    :return Network: Protein-protein interaction network with information on differential expression.
    """
    # Load and reify the graph
    bel_graph = bel_graph_loader(bel_graph_path)
    print('utils.generate_bel_network Test collapsing network before reifying')
    collapse_to_genes(bel_graph)
    print(bel_graph is None)
    reified = reify_bel_graph(bel_graph)

    if disease_associations_path is not None:
        if current_disease_ids_path:
            current_disease_ids = parse_disease_ids(current_disease_ids_path)
        else:
            current_disease_ids = set()
        disease_associations = parse_disease_associations(disease_associations_path,
                                                          current_disease_ids)
    else:
        disease_associations = None

    # Build an undirected weighted graph with the remaining interactions based on Entrez gene IDs
    network = NetworkNx(
        reified.to_undirected(),
        max_adj_p=max_adj_p,
        max_l2fc=max_log2_fold_change,
        min_l2fc=min_log2_fold_change,
    )
    network.set_up_network(dge_list, disease_associations=disease_associations)

    return network


def write_gat2vec_input_files(
        network: NetworkNx,
        targets: List[str],
        home_dir: str,
        assoc_score: Optional[Dict] = None
):
    """Write the input files for gat2vec tool.

    :param network: Network object with attributes overlayed on it.
    :param targets:
    :param home_dir:
    :param assoc_score:
    """
    network.write_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "graph"))

    attribute_network = AttributeFileGenerator(network)
    attribute_network.write_attribute_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "na"))

    labeled_network = LabeledNetworkGenerator(network)
    labeled_network.write_index_labels(
        targets,
        gat2vec_paths.get_labels_path(home_dir),
        sample_scores=assoc_score
    )


def add_disease_attribute(graph: nx.Graph, att_mappings: Dict):
    """Add the phenotypes to the Base Entities as attributes."""
    for node in graph:
        if isinstance(node, CentralDogma) and node.name in att_mappings:
            graph.nodes[node]['phenotypes'] = [phtype for _, phtype in att_mappings[node.name]]


def add_dge_attribute(graph: nx.Graph, gene_dict: Dict):
    """ """
    max_adjp = 0.05
    min_l2fc, max_l2fc = -0.5, 0.5
    dge = {
        g.entrez_id or g.symbol: g.log2_fold_change
        for g
        in gene_dict
        if g.padj < max_adjp
    }

    for node in graph:
        if isinstance(node, CentralDogma) and node['name'] in dge.keys():
            if dge[node['name']] > max_l2fc:
                graph.nodes[node]['deg'] = True
                graph.nodes[node]['upreg'] = True
            elif dge[node['name']] < min_l2fc:
                graph.nodes[node]['deg'] = True
                graph.nodes[node]['downreg'] = True


def get_significantly_differentiated(gene_list: List[Gene], max_adjp: float):
    """Returns a dictionary only with significantly differentially expressed genes from the gene list."""
    max_adjp = max_adjp or 0.05

    dge = {g.entrez_id or g.symbol: g.log2_fold_change
           for g in gene_list
           if g.padj < max_adjp}

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


def parse_gene_list(path: str, network: NetworkNx) -> list:
    """Parse a list of genes and return them if they are in the network.

    :param str path: The path of input file.
    :param Graph network: The graph with genes as nodes.
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
        network.graph.nodes[node]['name']
        for node
        in network.graph
        if 'name' in network.graph.nodes[node]
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


# TODO delete this?
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


def bel_graph_loader(from_dir: str) -> BELGraph:
    """Obtains a combined BELGraph from all the BEL documents in one folder.

    :param from_dir: The folder with the BEL documents.
    :return: A corresponding BEL Graph.
    """
    logger.info("Loading BEL Graph.")
    files = [
        join(from_dir, file)
        for file
        in listdir(from_dir)
        if isfile(join(from_dir, file))
    ]
    bel_files = [file for file in files if file[-4:].lower() == '.bel']
    bel_graphs = [from_path(file) for file in bel_files]
    return union(bel_graphs)


class AttributeFileGenerator:
    """Mimic encapsulation of a bipartite attribute network for Gat2Vec from a
    reified OpenBEL network."""

    def __init__(self, network: NetworkNx):
        """Initializes the network object.

        :param graph: The graph with nodes already converted to integers.
        """
        self.graph = network.graph
        self.network = network

    def write_attribute_adj_list(self, outpath: str):
        """Write the bipartite attribute graph to a file.

        :param str outpath: Path to the output file.
        """
        att_mappings = self.get_attribute_mappings()

        with open(outpath, mode="w") as file:
            for k, v in att_mappings.items():
                print("{} {}".format(k, " ".join(str(e) for e in v)), file=file)

    def get_attribute_mappings(self) -> Dict:
        att_ind_start = self.network.graph.number_of_nodes()

        att_mappings = defaultdict(list)
        att_ind_end = self._add_differential_expression_attributes(att_ind_start, att_mappings)
        att_ind_end = self._add_predicate_attributes(att_ind_end, att_mappings)
        if [node for node in self.network.graph if "associated_diseases" in self.network.graph.nodes[node]]:
            self._add_disease_association_attributes(att_ind_end, att_mappings)

        return att_mappings

    def _add_predicate_attributes(self, att_ind_start, att_mappings):
        """Add predicate information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :return: New index position.
        """
        predicate_mappings = self.get_predicate_mappings(att_ind_start)

        for node in self.network.graph:
            if (
                    type(self.network.graph.nodes[node]['bel_node']) == int and
                    self.network.graph.nodes[node]['label'] in predicate_mappings
            ):
                test = predicate_mappings[self.network.graph.nodes[node]['label']]
                att_mappings[node].append(test)
        return att_ind_start + len(predicate_mappings)

    def get_predicate_mappings(self, att_ind_start):
        """Get a dictionary of enumerations for predicates.

        :param int att_ind_start: Starting index for enumeration.
        :return: Dictionary of predicate, number pairs.
        """
        all_predicate_labels = self.get_all_unique_predicates()

        predicate_enum = enumerate(all_predicate_labels, start=att_ind_start)
        predicate_mappings = {}
        for num, pre in predicate_enum:
            predicate_mappings[pre] = num
        return predicate_mappings

    def get_all_unique_predicates(self):
        """Get all unique predicates that are known to the network.

        :return: All unique predicate identifiers.
        """
        # get unique elements
        all_predicate_labels = {
            self.network.graph.nodes[node]['label']
            for node
            in self.network.graph
            if type(self.network.graph.nodes[node]['bel_node']) == int
        }
        # remove None values and return as list
        return [lst for lst in all_predicate_labels if lst]

    def _add_disease_association_attributes(self, att_ind_start, att_mappings):
        """Add disease association information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        """
        disease_mappings = self.get_disease_mappings(att_ind_start)
        #
        for node in self.network.graph:
            if ('associated_diseases' in self.network.graph.nodes[node] and
                    self.network.graph.nodes[node]['associated_diseases']):
                assoc_disease_ids = [
                    disease_mappings[disease]
                    for disease
                    in self.network.graph.nodes[node]['associated_diseases']
                ]
                att_mappings[node].extend(assoc_disease_ids)

    def get_disease_mappings(self, att_ind_start):
        """Get a dictionary of enumerations for diseases.

        :param int att_ind_start: Starting index for enumeration.
        :return: Dictionary of disease, number pairs.
        """
        all_disease_ids = self.get_all_unique_diseases()

        disease_enum = enumerate(all_disease_ids, start=att_ind_start)
        disease_mappings = {}
        for num, dis in disease_enum:
            disease_mappings[dis] = num
        return disease_mappings

    def get_all_unique_diseases(self):
        """Get all unique diseases that are known to the network.

        :return: All unique disease identifiers.
        """
        all_disease_ids = [
            self.network.graph.nodes[node]['associated_diseases']
            for node
            in self.network.graph
            if ('associated_diseases' in self.network.graph.nodes[node] and
                self.network.graph.nodes[node]['associated_diseases'])
        ]
        # remove None values from list
        all_disease_ids = [lst for lst in all_disease_ids if lst is not None]
        # flatten list of lists, get unique elements
        all_disease_ids = set([
            phe
            for _list in all_disease_ids
            for phe in _list
            if phe
        ])

        return list(all_disease_ids)

    def _add_differential_expression_attributes(self, att_ind_start, att_mappings):
        """Add differential expression information to the attribute mapping dictionary.

        :param int att_ind_start: Start index for enumerating the attributes.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :return: End index for attribute enumeration.
        """
        up_regulated_ind = self.network.get_differentially_expressed_genes('up_regulated')
        down_regulated_ind = self.network.get_differentially_expressed_genes('down_regulated')
        rest_ind = self.network.get_differentially_expressed_genes('not_diff_expressed')

        self._add_attribute_values(att_ind_start, att_mappings, rest_ind)
        self._add_attribute_values(att_ind_start + 1, att_mappings, up_regulated_ind)
        self._add_attribute_values(att_ind_start + 2, att_mappings, down_regulated_ind)
        return att_ind_start + 3

    @staticmethod
    def _add_attribute_values(value, att_mappings, indices):
        """Add an attribute value to the given vertices.

        :param int value: Attribute value.
        :param dict att_mappings: Dictionary of mappings between vertices and enumerated attributes.
        :param list indices: Indices of the vertices.
        """
        for i in indices:
            att_mappings[i].append(value)


class LabeledNetworkGenerator:

    def __init__(self, network: NetworkNx):
        self.network = network

    def write_index_labels(self, targets, output_path, sample_scores: dict = None):
        """Write the mappings between vertex indices and labels(target vs. not) to a file.

        :param list targets: List of known targets.
        :param str output_path: Path to the output file.
        :param str sample_scores: Sample scores from OpenTarget.
        """
        label_mappings = self.get_index_labels(targets)

        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                if sample_scores:
                    if self.network.graph.nodes[k]["name"] in sample_scores:
                        score = self._convert_score_to_weight(
                            v,
                            sample_scores[self.network.graph.nodes[k]["name"]]
                        )
                        print(k, v, score, sep='\t', file=file)
                    else:
                        print(k, v, '1.', sep='\t', file=file)
                else:
                    print(k, v, sep='\t', file=file)

    def get_index_labels(self, targets):
        """Get the labels(known target/not) mapped to indices.

        :param targets: List of known targets
        :return: Dictionary of index-label mappings
        """
        gene_ind = {
            self.network.graph.nodes[node]['name']: node
            for node
            in self.network.graph
            if ('name' in self.network.graph.nodes[node])
        }
        target_ind = [node for name, node in gene_ind.items() if name in targets]
        rest_ind = [node for name, node in gene_ind.items() if name not in targets]
        label_mappings = {i: 1 for i in target_ind}
        label_mappings.update({i: 0 for i in rest_ind})
        return label_mappings

    @staticmethod
    def _convert_score_to_weight(label: int, score: float) -> float:
        """Convert the association score into a weight for the weighted classification. If the
        label is positive the weight is the score. If negative, the weight is 1 - score. This means,
        a high score for a negative label will imply some uncertainty about it being a target.

        :param label: 1 for positive, 0 for negative.
        :param score: The association score.
        :return: The weight.
        """
        if label:
            return score
        else:
            return 1 - score
