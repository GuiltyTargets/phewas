# -*- coding: utf-8 -*-

"""This module creates the networks required to predict links."""

from collections import defaultdict
import logging
import re
import sys
from typing import Any, Generator, List, Optional, Set

from igraph import EdgeSeq, Graph, Vertex
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from guiltytargets.gat2vec import gat2vec_paths
from guiltytargets.ppi_network_annotation.parsers import parse_ppi_graph
from guiltytargets.ppi_network_annotation.model import AttributeNetwork, Network
from guiltytargets_phewas.constants import lfc_cutoff
from guiltytargets_phewas.parsers import parse_disease_drug_graph, parse_disease_gene_graph, parse_gene_drug_graph

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
fh = logging.FileHandler('hetero_network.log')
fh.setLevel(logging.DEBUG)
logger.addHandler(ch)
logger.addHandler(fh)

drugbank_regex = re.compile('^DB\d{5}$')
disease_regex = re.compile('(OMIM|MESH|DOID|UMLS_CUI):')
# entrez_regex = re.compile('^\d+$')  <- Some genes are represented as Ensembl code or gene symbols.
# ENSP, ENSG, gene_symbol or entrez


def merge_by_attribute(g1: Graph, g2: Graph, v_attr_name='name') -> Graph:
    """Merge two graphs by one of the attributes. Copyied from https://github.com/igraph/igraph/issues/876.

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
    edge_list_by_attribute_2 = set([tuple(sorted([v_attr_2[u], v_attr_2[v]])) for u, v in g2.get_edgelist()])

    edge_list_by_attribute_merged = edge_list_by_attribute_1.union(edge_list_by_attribute_2)

    v_attr_merged = sorted(list(set(g2.vs[v_attr_name]).union(set(g1.vs[v_attr_name]))))

    attribute_to_ind = {v_attr_merged: i for i, v_attr_merged in enumerate(v_attr_merged)}
    edge_list_merged = [(attribute_to_ind[i], attribute_to_ind[j]) for i, j in edge_list_by_attribute_merged]

    graph_merged = Graph(edge_list_merged)
    graph_merged.vs[v_attr_name] = v_attr_merged

    # Include attributes that are in both g_1 and g_2. If different attribute values are present in a vertex,
    # then the one of g_1 is used
    for attr_name_other in set(g2.vs.attributes()).intersection(set(g1.vs.attributes())).difference([v_attr_name]):
        attr_other = dict()
        for v in g2.vs():
            attr_other[attribute_to_ind[v[v_attr_name]]] = v[attr_name_other]

        for v in g1.vs():
            attr_other[attribute_to_ind[v[v_attr_name]]] = v[attr_name_other]
        graph_merged.vs[attr_name_other] = [attr_other[i] for i in range(graph_merged.vcount())]

    return graph_merged


class HeterogeneousNetwork(Network):

    def __init__(
            self,
            ppi_graph: Graph,
            disease_gene_graph: Graph,
            disease_drug_graph: Graph,
            gene_drug_graph: Graph,
            infer_disease_target_assoc: bool = False
    ):
        """"""
        max_l2fc, min_l2fc = -1 * lfc_cutoff, lfc_cutoff
        super().__init__(ppi_graph, max_l2fc=max_l2fc, min_l2fc=min_l2fc)
        # count the number of merges for debug reasons
        if logger.getEffectiveLevel() == logging.DEBUG:
            drug1 = set(disease_drug_graph.vs.select(self._match_drug)['name'])
            drug2 = set(gene_drug_graph.vs.select(self._match_drug)['name'])
            merged = drug1.intersection(drug2)
            logger.debug(f'Merged drugs: {len(merged)}. Disease network: {len(drug1)}. Gene network: {len(drug2)}.')
            genes1 = set(disease_gene_graph.vs.select(self._match_gene)['name'])
            genes2 = set(gene_drug_graph.vs.select(self._match_gene)['name'])
            merged = genes1.intersection(genes2)
            logger.debug(f'Merged genes: {len(merged)}. Disease network: {len(genes1)}. Drug network: {len(genes2)}.')
            diseases1 = set(disease_gene_graph.vs.select(self._match_disease)['name'])
            diseases2 = set(disease_drug_graph.vs.select(self._match_disease)['name'])
            doid1 = [x for x in diseases1 if x.startswith('DOID:')]
            logger.debug(f'DOIDs in disease gene: {len(doid1)}')
            doid2 = [x for x in diseases2 if x.startswith('DOID:')]
            logger.debug(f'DOIDs in disease drug: {len(doid2)}')
            merged = diseases1.intersection(diseases2)
            logger.debug(f'Merged diseases: {len(merged)}. Gene network: {len(diseases1)}. '
                         f'Drug network: {len(diseases2)}.')

        self.graph = merge_by_attribute(ppi_graph, gene_drug_graph)
        self.graph = merge_by_attribute(self.graph, disease_drug_graph)
        self.graph = merge_by_attribute(self.graph, disease_gene_graph)
        self.set_node_types()
        if infer_disease_target_assoc:
            self._add_transitive_relations()

    def set_node_types(self):
        """Assign the types of nodes based on their ids."""
        self.graph.vs['type'] = 'gene'  # What is not disease or a drug, is a gene.
        self.graph.vs.select(self._match_disease)['type'] = 'disease'
        self.graph.vs.select(self._match_drug)['type'] = 'drug'

    def write_gat2vec_input_files(
            self,
            home_dir: str,
            disease_id: str,
            use_dge_data: bool = True,
            filter_pleiotropic_targets: bool = False
    ):
        """"""
        self.write_adj_list(
            path=gat2vec_paths.get_adjlist_path(home_dir, "graph"),
            ignore_disease=disease_id
        )
        attrib_network = self.get_attribute_network(use_dge_data)
        attrib_network.write_attribute_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "na"))
        labeled_network = self.get_labeled_network(filter_pleiotropic_targets)
        labeled_network.write_index_labels(disease_id, gat2vec_paths.get_labels_path(home_dir))

    def write_gat2vec_cv_split(
            self,
            home_dir: str,
            disease_id: str,
            use_dge_data: bool = True,
            filter_pleiotropic_targets: bool = True,
            n_splits: int = 5
    ) -> Generator[Any, None, None]:
        """"""
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True)

        attribute_network = self.get_attribute_network(use_dge_data)
        attribute_network.write_attribute_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "na"))
        labeled_network = self.get_labeled_network(filter_pleiotropic_targets)
        labels_dict = labeled_network.get_targets_index_labels(disease_id)
        logger.debug(f'Total vertices: {len(self.graph.vs)}')
        logger.debug(f'Gene vertices: {len(self.graph.vs.select(HeterogeneousNetwork.is_gene))}')
        logger.debug(f'True: {sum(1 for x in labels_dict.values() if x)}')
        logger.debug(f'False: {sum(1 for x in labels_dict.values() if x==0)}')

        df = pd.DataFrame(
            data=[
                (gene_vertex, label)
                for gene_vertex, label
                in labels_dict.items()
            ],
            columns=['vertex_id', 'labels']
        )

        for _, test_idx in skf.split(df.loc[:, 'vertex_id'], df.loc[:, 'labels']):
            logger.debug(f'Test size: {len(test_idx)}')
            cv_labels_dict = {
                row.vertex_id: row.labels for
                row in
                df.iloc[test_idx].itertuples(index=False)
            }
            targets_for_prediction = self.graph.vs.select(cv_labels_dict.keys())['name']
            logger.debug(f'Number of targets for prediction: {len(targets_for_prediction)}')
            logger.debug(f'Sample:')
            logger.debug(targets_for_prediction[:30])
            self.write_adj_list(
                path=gat2vec_paths.get_adjlist_path(home_dir, "graph"),
                ignore_targets=targets_for_prediction,
                ignore_disease=disease_id,
            )
            yield cv_labels_dict

    def _set_default_vertex_attributes(self) -> None:
        """Assign default values on attributes to all vertices."""
        genes = self.graph.vs.select(self.is_gene).indices
        logger.debug(f'# of genes: {len(genes)}')
        self.graph.vs["diff_expressed"] = None
        self.graph.vs["up_regulated"] = None
        self.graph.vs["down_regulated"] = None
        self.graph.vs(genes)["l2fc"] = 0
        self.graph.vs(genes)["padj"] = 0.5
        self.graph.vs(genes)["symbol"] = ''
        self.graph.vs(genes)["diff_expressed"] = False
        self.graph.vs(genes)["up_regulated"] = False
        self.graph.vs(genes)["down_regulated"] = False

    @staticmethod
    def is_drug(v: Vertex) -> bool:
        """Tests if a vertex is a drug node. """
        return v.attributes()['type'] == 'drug'

    @staticmethod
    def is_gene(v: Vertex) -> bool:
        """Tests if a vertex is a gene node. """
        return v.attributes()['type'] == 'gene'

    @staticmethod
    def is_disease(v: Vertex) -> bool:
        """Tests if a vertex is a disease node. """
        return v.attributes()['type'] == 'disease'

    @staticmethod
    def _match_drug(v: Vertex) -> bool:
        """Tests if a vertex name matches the regular expression for Drugbank ids. """
        return bool(drugbank_regex.match(v.attributes()["name"]))

    @staticmethod
    def _match_gene(v: Vertex) -> bool:
        """Tests if a vertex name matches doesn't match any other regex. """
        return not (HeterogeneousNetwork._match_disease(v) or HeterogeneousNetwork._match_drug(v))

    @staticmethod
    def _match_disease(v: Vertex) -> bool:
        """Tests if a vertex name matches the regular expression for disease ids."""
        try:
            disease_regex.match(v.attributes()["name"])
        except TypeError:
            print(v.attributes()["name"])
            print(type(v.attributes()["name"]))
            sys.exit()
        return bool(disease_regex.match(v.attributes()["name"]))

    def _is_significantly_differentiated(self, v: Vertex) -> bool:
        """In a heterogeneous network, the vertex has to be tested if it is a gene before testing if it's up or
        down regulated."""
        return self.is_gene(v) and v.attributes()['padj'] < self.max_adj_p

    def get_disease_gene_assoc(self, from_disease: str = None) -> EdgeSeq:
        """All the edges between a disease and a gene."""
        genes = self.graph.vs.select(self.is_gene).indices
        if from_disease:
            disease = [self.get_index_for_disease(from_disease)]
        else:
            disease = self.graph.vs.select(self.is_disease).indices
        return self.graph.es.select(_between=(genes, disease))

    def write_adj_list(
            self,
            path: str,
            ignore_targets: Optional[List[str]] = None,
            ignore_disease: str = None,
    ) -> None:
        """Saves the graph as an adjacency list file. May choose to ignore edges from a disease to all targets, from
        a disease to a list of targets or from all diseases to

        :param path: Path where the file will be writen.
        :param ignore_targets: A list of targets to be ignored. Ignored if no disease is specified.
        :param ignore_disease: Specify a disease name to have all edges between its respective node and the
         targets ignored. Use None if no disease is to be ignored.
        :return:
        """
        if ignore_targets and ignore_disease:
            full_graph = self.get_full_graph()
            disease_idx = full_graph.vs.find(name_eq=ignore_disease).index
            genes = full_graph.vs.select(name_in=ignore_targets).indices
            to_remove = self.graph.es.select(_between=(genes, [disease_idx]))
            full_graph.delete_edges(to_remove)
            assert disease_idx == full_graph.vs.find(name_eq=ignore_disease).index, 'The index of the disease changed!'
            adj_list = full_graph.get_adjlist()
        else:
            adj_list = self.get_full_graph(ignore_disease).get_adjlist()
        with open(path, mode="w") as file:
            for i, line in enumerate(adj_list):
                print(i, *line, file=file)
        """import os
from guiltytargets_phewas.target_repositioning  import generate_heterogeneous_network
base = r'/home/bit/lacerda/data/'
ppi = os.path.join(base, 'STRING', 'string_entrez.edgelist')
dg = os.path.join(base, 'SNAP', 'DG-AssocMiner_miner-disease-gene.tsv.gz')
dc = os.path.join(base, 'SNAP', 'DCh-Miner_miner-disease-chemical.tsv.gz')
gc = os.path.join(base, 'SNAP', 'ChG-Miner_miner-chem-gene.tsv.gz')
net = generate_heterogeneous_network(ppi, dg, dc, gc, True)
net.write_adj_list('sai.txt', ignore_disease='DOID:2377', ignore_targets=['3575'])
net.write_adj_list('sai.txt', ignore_disease='DOID:2377')
len(list(net.get_full_graph('DOID:2377').es))
with open('sai.txt') as f:
   print(len(f.read().splitlines()))


# a7806244
# b7805077
"""

    def get_full_graph(self, ignore_disease_gene: str = None) -> Graph:
        """Creates a copy of the graph. If a disease is specified, ignore the edges between all targets to this disease.

        :param ignore_disease_gene: The description of a disease to be ignored.
        :return: A copy of the graph.
        """
        if ignore_disease_gene:
            graph_copy = self.graph.copy()
            ignored_edges = self.get_disease_gene_assoc(ignore_disease_gene)
            graph_copy.delete_edges(ignored_edges)
            return graph_copy
        else:
            return self.graph.copy()

    def _write_drug_genes_list(self, path: str):
        """Creates the disease gene association file as an edgelist file."""
        adj_dict = defaultdict(list)

        # Create the adjacency list between disease and gene vertices.
        for e in self.get_disease_gene_assoc():
            v1, v2 = e.tuple
            adj_dict[v1].append(v2)
            adj_dict[v2].append(v1)

        with open(path, mode='w') as file:
            for i, line in adj_dict.items():
                print(i, *line, file=file)

    def _write_disease_genes_association(self, path: str):
        """Creates a file with connected vertex pairs."""
        with open(path, mode='w') as file:
            for e in self.get_disease_gene_assoc():
                print(" ".join(str(x) for x in e.tuple), file=file)

    def get_labeled_network(self, filter_pleiotropic_targets: bool = False):
        """Generates labels according to a single disease."""
        return HeterogeneousLabelGenerator(self, filter_pleiotropic_targets=filter_pleiotropic_targets)

    def get_attribute_network(self, use_dge_data: bool = True):
        """Generates labels according to a single disease."""
        if use_dge_data:
            return AttributeNetwork(self)
        else:
            return HeterogenousAttributeGenerator(self)

    def find_pleiotropic_genes(self) -> Set[int]:
        """Finds genes that are associated with more than one disease."""
        genes = self.graph.vs.select(self.is_gene).indices
        diseases = self.graph.vs.select(self.is_disease).indices
        gene_neighbors = self.graph.neighborhood(genes)
        neighbor_diseases = [
            set(n_list).intersection(diseases)
            for n_list
            in gene_neighbors
        ]
        return set(
            gene
            for gene, assoc_diseases
            in zip(genes, neighbor_diseases)
            if len(assoc_diseases) > 1
        )

    def _add_transitive_relations(self):
        """This checks unchecked disease-gene relations inferred by disease-drug and drug-gene relations (If drug1 is
         indicated to disease1 and drug1 targets gene1, gene1 is associated with disease1)."""
        graph_copy = self.graph.copy()

        # Remove the edges connecting disease and genes...
        d_list = graph_copy.vs.select(self.is_disease).indices
        g_list = graph_copy.vs.select(self.is_gene).indices
        dg_edges = graph_copy.es.select(_between=(d_list, g_list))
        graph_copy.delete_edges(dg_edges)
        # ... so all disease-gene paths with order 2 have a disease-drug-gene metapath, the ones we want to infer.
        neighborhood = graph_copy.neighborhood(d_list, order=2)
        filtered_neighborhood = [
            list(set(n_list).intersection(g_list))
            for n_list
            in neighborhood
        ]
        new_e_list = [
            (d_idx, neighbor_gene)
            for d_idx, n_list
            in zip(d_list, filtered_neighborhood)
            for neighbor_gene in n_list
        ]
        self.graph.add_edges(new_e_list)

    def get_index_for_disease(self, disease):
        """"""
        return self.graph.vs.find(name_eq=disease).index


class HeterogenousAttributeGenerator(AttributeNetwork):
    """"""

    def get_attribute_mappings(self):
        """Get a dictionary of mappings between vertices and enumerated attributes. This considers the node types as
        attributes.

        :return: Dictionary of mappings between vertices and enumerated attributes.
        """
        att_ind_start = len(self.graph.vs)
        att_mappings = {
            x: [att_ind_start + 1]
            for x
            in self.graph.vs.indices
        }
        return att_mappings


class HeterogeneousLabelGenerator:
    """This labeled network will write pairs of connected nodes that will be ranked and the pairs that will
    be disconnected."""

    def __init__(self, network: HeterogeneousNetwork, filter_pleiotropic_targets: bool = False):
        """Initialize the network object.

        :param network: A heterogeneous network (containing targets, diseases and drugs) annotated with differential
        gene expression and disease association.
        """
        self.network = network
        self._filter_pleiotropic_targets = filter_pleiotropic_targets

    def write_index_labels(self, disease: str, output_path):
        """Write the mappings between the specified disease vertex to the genes (connected/ disconnected) to a file.

        :param disease: Disease being assessed.
        :param output_path: Path to the output file.
        :return:
        """
        disease_idx = self.network.graph.vs.find(name_eq=disease).index

        label_mappings = self.get_targets_index_labels(disease)
        
        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                print(disease_idx, k, v, sep='\t', file=file)

    def get_targets_index_labels(self, from_disease):
        """Get the labels(connected/not connected) mapped to indices. """
        disease_idx = self.network.get_index_for_disease(from_disease)
        # disease_idx = self.network.graph.vs.find(name_eq=from_disease).index

        neighbors_ind = self.network.graph.neighbors(disease_idx)
        neighbors_names = self.network.graph.vs(neighbors_ind)['name']

        genes = self.network.graph.vs.select(HeterogeneousNetwork.is_gene)

        connected_ind = self.network.graph.vs.select(HeterogeneousNetwork.is_gene, name_in=neighbors_names).indices
        if self._filter_pleiotropic_targets:
            connected_ind = self.network.find_pleiotropic_genes().intersection(connected_ind)

        label_mappings = {
            i: 1 if i in connected_ind else 0
            for i
            in genes.indices
        }
        return label_mappings


def generate_heterogeneous_network(
        ppi_graph_path,
        disease_gene_path: str,
        disease_drug_path: str,
        gene_drug_path: str,
        infer_disease_target_assoc: bool = True
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
        gene_drug_interactions,
        infer_disease_target_assoc=infer_disease_target_assoc
    )
    return network
