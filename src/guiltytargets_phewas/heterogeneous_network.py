# -*- coding: utf-8 -*-

"""This module predict links."""

import logging
import re
from typing import List, Optional

from igraph import EdgeSeq, Graph, Vertex

from guiltytargets.ppi_network_annotation.parsers import parse_ppi_graph
from guiltytargets.ppi_network_annotation.model import Gene
from .parsers import parse_disease_drug_graph, parse_disease_gene_graph, parse_gene_drug_graph

from .utils import uniprot_id_to_entrez_id_converter

logger = logging.getLogger(__name__)

drugbank_regex = re.compile('^DB\d{5}$')
disease_regex = re.compile('(OMIM|MESH):')
entrez_regex = re.compile('^\d+$')
"""uniprot_regex = re.compile('^([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})' + '|' +
                             '([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$')"""


def merge_by_attribute(g1: Graph, g2: Graph, v_attr_name='name') -> Graph:
    """Merge two graphs by one of the attributes. Ignore other existing attributes.

    :param g1: Graph 1.
    :param g2: Graph 2.
    :param v_attr_name: The name of the attribute to be used to merge.
    :return: The merged graph.
    """
    assert len(set(g1.vs[v_attr_name])) == len(g1.vs[v_attr_name]), "Merging attribute must be unique"
    assert len(set(g2.vs[v_attr_name])) == len(g2.vs[v_attr_name]), "Merging attribute must be unique"

    v_attr_1 = g1.vs[v_attr_name]
    v_attr_2 = g2.vs[v_attr_name]

    edge_list_by_attribute_1 = set([tuple(sorted([v_attr_1[u], v_attr_1[v]])) for u, v in g1.get_edgelist()])
    edge_list_by_attribute_2 = set([tuple(sorted([v_attr_1[u], v_attr_1[v]])) for u, v in g2.get_edgelist()])

    edge_list_by_attribute_merged = edge_list_by_attribute_1.union(edge_list_by_attribute_2)

    v_attr_merged = sorted(list(set(v_attr_2).union(set(v_attr_1))))

    attribute_to_ind = {v_attr_merged: i for i, v_attr_merged in enumerate(v_attr_merged)}
    edge_list_merged = [(attribute_to_ind[i], attribute_to_ind[j]) for i, j in edge_list_by_attribute_merged]

    graph_merged = Graph(edge_list_merged)
    graph_merged.vs[v_attr_name] = v_attr_merged

    return graph_merged


class HeterogeneousNetwork:

    gene = 'entrez'
    disease = 'disease'
    drug = 'drug'
    max_adj_p = 0.5
    min_l2fc = 1.5
    max_l2fc = -1.5

    def __init__(
            self,
            ppi_graph: Graph,
            disease_gene_graph: Graph,
            disease_drug_graph: Graph,
            gene_drug_graph: Graph,
    ):
        """"""
        # count the number of merges for debug reasons
        genes1 = set(disease_gene_graph.vs.select(self._is_gene)['name'])
        genes2 = set(gene_drug_graph.vs.select(self._is_gene)['name'])
        merged = genes1.intersection(genes2)
        logger.info(f'Merged genes: {len(merged)}. Disease network: {len(genes1)}. Drug network: {len(genes2)}.')
        diseases1 = set(disease_gene_graph.vs.select(self._is_disease)['name'])
        diseases2 = set(disease_drug_graph.vs.select(self._is_disease)['name'])
        merged = diseases1.intersection(diseases2)
        logger.info(f'Merged genes: {len(merged)}. Gene network: {len(diseases1)}. Drug network: {len(diseases2)}.')
        self.graph = merge_by_attribute(gene_drug_graph, disease_gene_graph)
        self.graph = merge_by_attribute(self.graph, disease_drug_graph)
        self.graph = merge_by_attribute(self.graph, ppi_graph)

    def setup_network(
            self,
            genes: Optional[List[Gene]] = None,
    ):
        """"""
        self._add_vertex_attributes(genes)
        self._add_node_type_information()

    def _set_default_vertex_attributes(self):
        """Assign default values on attributes to all vertices."""
        self.graph.vs["l2fc"] = 0
        self.graph.vs["padj"] = 0.5
        self.graph.vs["diff_expressed"] = False
        self.graph.vs["up_regulated"] = False
        self.graph.vs["down_regulated"] = False

    def _add_vertex_attributes(self, genes: Optional[List[Gene]] = None):
        self._set_default_vertex_attributes()

    def _add_node_type_information(self):
        """"""
        self.graph.vs.select(self._is_drug)['type'] = 'Drug'
        self.graph.vs.select(self._is_disease)['type'] = 'Disease'
        self.graph.vs.select(self._is_gene)['type'] = 'Gene'

    @staticmethod
    def _is_drug(v: Vertex) -> bool:
        """Tests if a vertex is a drug node. """
        return bool(drugbank_regex.match(v.attributes()["name"]))

    @staticmethod
    def _is_gene(v: Vertex) -> bool:
        """Tests if a vertex is a gene node. """
        return bool(entrez_regex.match(v.attributes()["name"]))

    @staticmethod
    def _is_disease(v: Vertex) -> bool:
        """Tests if a vertex is a disease node. """
        return bool(disease_regex.match(v.attributes()["name"]))

    def get_drug_gene_assoc(self) -> EdgeSeq:
        """All the edges between a drug and a gene."""
        genes = self.graph.vs.select(self._is_gene).indices
        drugs = self.graph.vs.select(self._is_drug).indices
        return self.graph.es.select(_between=(genes, drugs))


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
