# -*- coding: utf-8 -*-

"""This module contains the extension to the class Network from PPI annotation framework."""

import logging
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
from pybel.dsl import gene, protein, rna

from guiltytargets.ppi_network_annotation import Gene

__all__ = [
    'NetworkNx',
]

logger = logging.getLogger(__name__)


class NetworkNx:
    """Encapsulate a PPI network with differential gene expression, phenotypes and disease association annotation."""

    def __init__(
            self,
            ppi_graph: nx.Graph,
            max_adj_p: Optional[float] = None,
            max_l2fc: Optional[float] = None,
            min_l2fc: Optional[float] = None,
    ) -> None:
        """Initialize the network object.

        :param ppi_graph: A graph of protein interactions.
        :param max_adj_p: Maximum value for adjusted p-value, used for calculating differential expression.
        :param max_l2fc: Maximum value for log2 fold change, used for calculating down regulation.
        :param min_l2fc: Minimum value for log2 fold change, used for calculating up regulation.
        """

        logger.info("Initializing Network")

        self.max_adj_p = max_adj_p or 0.05
        self.max_l2fc = max_l2fc or -1.0
        self.min_l2fc = min_l2fc or +1.0

        self.graph = ppi_graph.copy()

    def set_up_network(
            self,
            genes: List[Gene],
            gene_filter: bool = False,
            disease_associations: Optional[Dict] = None
    ) -> None:
        """Set up the network.

         Filter genes out if requested and add attributes to the vertices.

        :param genes: A list of Gene objects.
        :param gene_filter: Removes all genes that are not in list <genes> if True.
        :param disease_associations: Diseases associated with genes.
        """
        if gene_filter:
            self.filter_genes([gene.entrez_id for gene in genes])
        self._add_vertex_attributes(genes, disease_associations)

    def filter_genes(self, relevant_entrez: list) -> None:
        """Filter out the genes that are not in list relevant_entrez.

        :param list relevant_entrez: Entrez IDs of genes which are to be kept.
        """
        logger.info("In filter_genes()")
        raise Exception('Not ready to filter genes using NetworkX')

    def _add_vertex_attributes(
            self,
            genes: List[Gene],
            disease_associations: Optional[dict] = None,
    ) -> None:
        """Add attributes to vertices.

        :param genes: A list of genes containing attribute information.
        """
        self._set_default_vertex_attributes()
        self._add_vertex_attributes_by_genes(genes)

        # set the attributes for up-regulated and down-regulated genes
        up_regulated, down_regulated = 0, 0
        for node in self.graph:
            cur_node = self.graph.nodes[node]
            if isinstance(cur_node['bel_node'], (protein, rna, gene)):
                if cur_node['padj'] < self.max_adj_p:
                    if cur_node['l2fc'] > self.max_l2fc:
                        cur_node['diff_expressed'] = True
                        cur_node['up_regulated'] = True
                        up_regulated += 1
                    elif cur_node['l2fc'] > self.min_l2fc:
                        cur_node['diff_expressed'] = True
                        cur_node['down_regulated'] = True
                        down_regulated += 1

        # add disease associations
        self._add_disease_associations(disease_associations)

        logger.info("Number of all differentially expressed genes is: {}".
                    format(up_regulated + down_regulated))

    def _set_default_vertex_attributes(self):
        """Initializes each vertex with its required attributes."""
        for node in self.graph:
            if isinstance(self.graph.nodes[node]['bel_node'], (protein, rna, gene)):
                self.graph.nodes[node]["l2fc"] = 0
                self.graph.nodes[node]["padj"] = 0.5
                self.graph.nodes[node]["symbol"] = self.graph.nodes[node]['bel_node']["name"]
                self.graph.nodes[node]["diff_expressed"] = False
                self.graph.nodes[node]["up_regulated"] = False
                self.graph.nodes[node]["down_regulated"] = False

    def _add_vertex_attributes_by_genes(self, gene_list: List[Gene]) -> None:
        """Assign values to attributes on vertices.

        :param gene_list: A list of Gene objects from which values will be extracted.
        """
        gene_dict = {
            g.entrez_id: {
                'l2fc': g.log2_fold_change,
                'symbol': g.symbol,
                'padj': g.padj
            }
            for g
            in gene_list
        }
        for node in self.graph:
            cur_node = self.graph.nodes[node]
            if (isinstance(cur_node['bel_node'], (protein, rna, gene)) and
                    cur_node['bel_node']['name'] in gene_dict):
                cur_gene = gene_dict[cur_node['bel_node']['name']]
                cur_node['l2fc'] = cur_gene['l2fc']
                cur_node['symbol'] = cur_gene['symbol']
                cur_node['padj'] = cur_gene['padj']

    def _add_disease_associations(self, disease_associations: dict) -> None:
        """Add disease association annotation to the network.

        :param disease_associations: Dictionary of disease-gene associations.
        """
        if disease_associations is not None:
            for node in self.graph:
                cur_node = self.graph.nodes[node]
                if ('bel_node' in cur_node and
                        isinstance(cur_node['bel_node'], (protein, rna, gene)) and
                        cur_node['bel_node']['name'] in disease_associations):
                    cur_node['associated_diseases'] = disease_associations[cur_node['bel_node']['name']]

    def write_adj_list(self, path: str) -> None:
        """Write the network as an adjacency list to a file.

        :param path: File path to write the adjacency list.
        """
        nx.write_adjlist(self.graph, path)

    def get_attribute_from_indices(self, indices: list, attribute_name: str):
        """Get attribute values for the requested indices.

        :param indices: Indices of vertices for which the attribute values are requested.
        :param attribute_name: The name of the attribute.
        :return: A list of attribute values for the requested indices.
        """
        # TODO Update for NetworkX
        return list(np.array(self.graph.vs[attribute_name])[indices])

    def get_differentially_expressed_genes(self, diff_type: str) -> List:
        """Get up regulated, down regulated, all differentially expressed or not differentially expressed nodes.

        :param diff_type: One of `not_diff_expressed`, `diff_expressed`, `up_regulated`, `down_regulated`
        :return: A list of nodes corresponding to diff_type.
        """
        assert diff_type in [
            "not_diff_expressed", "diff_expressed", "up_regulated", "down_regulated"
        ], f"{diff_type} is not a valid type"
        nodes = []
        for node in self.graph:
            if {"diff_expressed", "up_regulated", "down_regulated"}.issubset(self.graph.nodes[node].keys()):
                if diff_type == "not_diff_expressed":
                    if not self.graph.nodes[node]["diff_expressed"]:
                        nodes.append(node)
                elif self.graph.nodes[node][diff_type]:
                    nodes.append(node)
        return nodes
