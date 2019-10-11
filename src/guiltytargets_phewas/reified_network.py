# -*- coding: utf-8 -*-

"""This module contains the extension to the class Network from PPI annotation framework."""

from collections import defaultdict
import logging
from typing import Dict, List, Optional

import networkx as nx
from pybel.dsl import CentralDogma

from guiltytargets.ppi_network_annotation import Gene

logger = logging.getLogger(__name__)

__all__ = [
    'NetworkNx',
    'LabeledNetworkGenerator',
    'AttributeFileGenerator',
]


class NetworkNx:
    """Encapsulate a PPI network with differential gene expression, phenotypes and disease association annotation."""

    def __init__(
            self,
            reified_graph: nx.Graph,
            max_adj_p: Optional[float] = None,
            max_l2fc: Optional[float] = None,
            min_l2fc: Optional[float] = None
    ) -> None:
        """Initialize the network object.

        :param reified_graph: A graph of protein interactions.
        :param max_adj_p: Maximum value for adjusted p-value, used for calculating differential expression.
        :param max_l2fc: Maximum value for log2 fold change, used for calculating down regulation.
        :param min_l2fc: Minimum value for log2 fold change, used for calculating up regulation.
        """

        logger.info("Initializing Network")

        self.max_adj_p = max_adj_p or 0.05
        self.max_l2fc = max_l2fc or -1.0
        self.min_l2fc = min_l2fc or +1.0

        self.graph = nx.relabel.convert_node_labels_to_integers(
            reified_graph,
            label_attribute='bel_node'
        )

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
            self.filter_genes([g.entrez_id for g in genes])
        self._add_vertex_attributes(genes, disease_associations)

    def filter_genes(self, relevant_entrez: list) -> None:
        """Filter out the genes that are not in list relevant_entrez.

        :param list relevant_entrez: Entrez IDs of genes which are to be kept.
        """
        logger.info("In filter_genes()")
        irrelevant = [node for node in self.graph if node not in relevant_entrez]
        self.graph.remove_nodes_from(irrelevant)
        # TODO Remove disconnected nodes
        raise Exception('Not ready to filter genes using NetworkX')

    def _add_vertex_attributes(
            self,
            genes: List[Gene],
            disease_associations: Optional[dict] = None
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
            if isinstance(cur_node['bel_node'], CentralDogma):
                if cur_node['padj'] < self.max_adj_p:
                    if cur_node['l2fc'] < self.max_l2fc:
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
            if isinstance(self.graph.nodes[node]['bel_node'], CentralDogma):
                self.graph.nodes[node]["l2fc"] = 0
                self.graph.nodes[node]["padj"] = 0.5
                self.graph.nodes[node]["name"] = self.graph.nodes[node]['bel_node']["name"]
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
            if (isinstance(cur_node['bel_node'], CentralDogma) and
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
                        isinstance(cur_node['bel_node'], CentralDogma) and
                        cur_node['bel_node']['name'] in disease_associations):
                    cur_node['associated_diseases'] = disease_associations[cur_node['bel_node']['name']]

    def write_adj_list(self, path: str) -> None:
        """Write the network as an adjacency list to a file.

        :param path: File path to write the adjacency list.
        """
        nx.write_adjlist(self.graph, path)

    def get_attribute_from_indices(self, indices: list, attribute_name: str) -> List:
        """Get attribute values for the requested indices.

        :param indices: Indices of vertices for which the attribute values are requested.
        :param attribute_name: The name of the attribute.
        :return: A list of attribute values for the requested indices.
        """
        # TODO Not sure how this works.
        gene_ind = {
            node: self.graph.nodes[node][attribute_name]
            for node
            in self.graph
            if attribute_name in self.graph.nodes[node]
        }

        return [
            attribute
            for node, attribute
            in gene_ind.items()
            if node in indices
        ]

    def get_differentially_expressed_genes(self, diff_type: str) -> List:
        """Get up regulated, down regulated, all differentially expressed or not
        differentially expressed nodes indices.

        :param diff_type: One of `not_diff_expressed`, `diff_expressed`, `up_regulated`, `down_regulated`
        :return: A list of nodes corresponding to diff_type.
        """
        assert diff_type in [
            "not_diff_expressed", "diff_expressed", "up_regulated", "down_regulated"
        ], f"{diff_type} is not a valid type"
        if diff_type == "not_diff_expressed":
            return [
                node
                for node
                in self.graph
                if ('diff_expressed' in self.graph.nodes[node] and
                    not self.graph.nodes[node]['diff_expressed'])
            ]
        else:
            return [
                node
                for node
                in self.graph
                if (diff_type in self.graph.nodes[node] and
                    self.graph.nodes[node][diff_type])
            ]

    def get_labeled_network(self):
        return LabeledNetworkGenerator(self)

    def get_attribute_network(self):
        return AttributeFileGenerator(self)

    def find_genes(self, gene_list: list, anno_type: str = "name"):
        gene_list = set(gene_list)
        symbols = [
            self.graph.nodes[node]['name']
            for node
            in self.graph
            if 'name' in self.graph.nodes[node]
        ]

        return list(gene_list.intersection(symbols))


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
                    isinstance(self.network.graph.nodes[node]['bel_node'], int) and
                    self.network.graph.nodes[node]['label'] in predicate_mappings
            ):
                att_mappings[node].append(predicate_mappings[self.network.graph.nodes[node]['label']])
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

    def get_all_unique_predicates(self) -> list:
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
        return list(filter(None, all_predicate_labels))

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

    def _add_go_terms_attributes(self):
        """"""
        # Get all unique GO terms
        #


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

    def write_files_for_g2v(
            self,
            targets: List[str],
            home_dir: str,
            assoc_score: Optional[Dict] = None
    ):
        ...

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
