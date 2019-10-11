# -*- coding: utf-8 -*-

"""Utility functions to run the algorithm."""

import logging
from os import listdir
from os.path import isfile, join
import time
import traceback
from typing import Dict, List, Optional, Set, Union
from urllib import parse, request

import bio2bel_phewascatalog
from igraph import Graph
import mygene
import networkx as nx
import pandas as pd
from opentargets import OpenTargetsClient
from pybel import BELGraph, from_path, union
from pybel.dsl import BaseEntity, CentralDogma
from pybel_tools.assembler.reified_graph.assembler import reify_bel_graph

from guiltytargets.gat2vec import gat2vec_paths
from guiltytargets.ppi_network_annotation.model.gene import Gene
from guiltytargets.ppi_network_annotation.parsers import (
    parse_disease_ids, parse_disease_associations
)
from .reified_network import AttributeFileGenerator, LabeledNetworkGenerator, NetworkNx

logger = logging.getLogger(__name__)


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
        collapse: str = None
):
    """Generate the protein-protein interaction network.

    :param bel_graph_path:
    :param dge_list:
    :param max_adj_p:
    :param max_log2_fold_change:
    :param min_log2_fold_change:
    :param current_disease_ids_path:
    :param disease_associations_path:
    :param collapse:
    :return: Interaction network with information on differential expression.
    """
    # Load and reify the graph
    bel_graph = bel_graph_loader(bel_graph_path)
    reified = reify_bel_graph(bel_graph, collapse=collapse)

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
        reified,
        max_adj_p=max_adj_p,
        max_l2fc=max_log2_fold_change,
        min_l2fc=min_log2_fold_change,
    )
    network.set_up_network(dge_list, disease_associations=disease_associations)

    return network


def write_gat2vec_input_files(
        network,
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


# TODO Renamed from download_for_disease. Moved to assembler. Delete.
def generate_targets_file(disease_id, outpath, anno_type: str = 'entrezgene') -> None:
    """Creates a disease list

    :param disease_id: EFO code from the disease.
    :param outpath:
    :param anno_type: `entrezgene` for Entrez Id or `symbol` for Gene symbol.
    :return:
    """
    ot = OpenTargetsClient()
    assoc = ot.get_associations_for_disease(
        disease_id,
        fields=['association_scoredatatypes', 'target.id']
    ).filter(
        datatype='known_drug'
    )
    ensembl_list = [a['target']['id'] for a in assoc]

    mg = mygene.MyGeneInfo()
    id_mappings = mg.getgenes(ensembl_list, fields=anno_type)

    with open(outpath, 'w+') as outfile:
        for mapping in id_mappings:
            if anno_type in mapping.keys():
                outfile.write(mapping[anno_type])
                outfile.write('\n')


# TODO Was renamed from get_association scores. Moved to input. Delete.
def generate_disease_gene_association_file(disease_id, outpath, anno_type: str = 'entrezgene'):
    """Obtain the association scores from the specified disease that are
    stored in the OpenTargets database.

    :param disease_id: The EFO code to the disease.
    :param outpath: The path to the file to be created.
    :param anno_type: `entrezgene` for Entrez Id or `symbol` for Gene symbol.
    :return:
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
    id_mappings = get_converter_to_entrez(ensembl_list)

    # Get the symbols and the scores
    ensembl_list = [
        (id_mappings[a['id']], a['score'])
        for a in assoc_simple
        if a['id'] in id_mappings
    ]

    with open(outpath, 'w+') as outfile:
        for symbol, score in ensembl_list:
            print(f'{symbol}\t{score}', file=outfile)


# TODO Moved to parsers. Delete.
def get_converter_to_entrez(
        query_list: List[str]
) -> Dict[str, str]:
    """Get a converter from a gene field to Entrez id.

    :param query_list: List of genes to query.
    :return: a dictionary with the query element as keys and entrez id as values.
    """
    mg = mygene.MyGeneInfo()
    gene_query = mg.querymany(query_list, scopes='symbol,entrezgene,ensembl.gene', species=9606, returnall=True)
    gene_query = [x for x in gene_query['out'] if 'entrezgene' in x]
    id_mappings = {
        g['query']: g['entrezgene']
        for g in gene_query
        if 'entrezgene' in g
    }

    logger.debug(f'get_converter_to_entrez - {len(query_list) - len(id_mappings)} not found.')
    logger.debug(f'get_converter_to_entrez - {len(id_mappings)} mappings')

    return id_mappings


# TODO moved this to parsers (GuiltyTargets). Delete.
def parse_gene_list_reified(path: str, graph: nx.Graph) -> list:
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
        node['name']
        for node
        in graph.nodes
        if isinstance(node, CentralDogma) and 'name' in node
    ]

    return list(genes.intersection(symbols))


# TODO moved this to parsers (GuiltyTargets). Delete.
def parse_gene_list(path: str, network) -> list:
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
    return network.find_genes(genes)


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


# TODO Moved to parsers. Delete.
def generate_phewas_file(out_path: str, anno_type: str = 'symbol'):
    """Create a file with gene-disease association from PheWAS data.

    :param out_path: Path where the file will be created.
    :param anno_type: `entrezgene` for Entrez Id or `symbol` for Gene symbol.
    :return:
    """
    phewas_manager = bio2bel_phewascatalog.Manager()
    pw_dict = phewas_manager.to_dict()
    if anno_type == 'entrezgene':
        to_entrez = get_converter_to_entrez(list(pw_dict.keys()))
        pw_dict = {
            entrez: pw_dict[symbol]
            for symbol, entrez
            in to_entrez.items()
        }
    with open(out_path, 'w') as file:
        for target, diseaselist in pw_dict.items():
            for disease in diseaselist:
                print(f'{target} {disease[1]}', file=file)


# TODO delete this only after notebook BEL2Gat2Vec(DGE) becomes obsolete
def get_phenotype_list(rbg: nx.Graph) -> List:
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


# TODO Moved to Parsers. Delete.
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


# TODO moved to converters. Delete.
def uniprot_id_to_entrez_id_converter(uniprot_id_list: Union[List[str], Set[str]]) -> Dict:
    """Converts a list of gene identifiers, from Uniprot to Entrez, using the Uniprot API.

    :param uniprot_id_list: The list of identifiers.
    :return: A converter dictionary from Uniprot id to Entrez id.
    """
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': 'ACC+ID',
        'to': 'P_ENTREZGENEID',
        'format': 'tab',
        'query': ' '.join(uniprot_id_list)
    }
    data = parse.urlencode(params)
    data = data.encode('utf-8')
    req = request.Request(url, data)
    with request.urlopen(req) as f:
        response = f.read()
    conv = {
        line.split('\t')[0]: line.split('\t')[1]
        for line
        in response.decode('utf-8').split('\n')
        if line
    }
    return conv


# to be used for link prediction
def merge_by_attribute(g1: Graph, g2: Graph, v_attr_name) -> Graph:
    """Merge two graphs by one of the attributes (`name` or `symbol`).

    :param g1: Graph 1.
    :param g2: Graph 2.
    :param v_attr_name: The name of the attribute to be used to merge.
    :return: The merged graph.
    """
    assert len(set(g1.vs[v_attr_name])) == len(g1.vs[v_attr_name]), "Merging attribute must be unique"
    assert len(set(g2.vs[v_attr_name])) == len(g2.vs[v_attr_name]), "Merging attribute must be unique"

    v_attr_1 = g1.vs[v_attr_name]
    v_attr_2 = g2.vs[v_attr_name]

    edge_list_by_attribute_1 = set([
        tuple(sorted([v_attr_1[u], v_attr_1[v]]) + [g1.es[(u, v)]['edge_type']])
        for u, v
        in g1.get_edgelist()
    ])
    edge_list_by_attribute_2 = set([
        tuple(sorted([v_attr_2[u], v_attr_2[v]]) + [g1.es[(u, v)]['edge_type']])
        for u, v
        in g2.get_edgelist()
    ])

    # Edges should not coincide (different edge types)
    edge_list_by_attribute_merged = edge_list_by_attribute_1.union(edge_list_by_attribute_2)

    v_attr_merged = sorted(list(set(v_attr_2).union(set(v_attr_1))))

    attribute_to_ind = {v_attr_merged: i for i, v_attr_merged in enumerate(v_attr_merged)}
    edge_list_merged = {
        (attribute_to_ind[i], attribute_to_ind[j]): e_type
        for i, j, e_type
        in edge_list_by_attribute_merged
    }

    graph_merged = Graph(list(edge_list_merged.keys()))
    for edge, e_type in edge_list_merged.items():
        graph_merged.es[edge]['edge_type'] = e_type
    graph_merged.vs[v_attr_name] = v_attr_merged

    # Include attributes that are in both g1 and g2. If different attribute values are present in a vertex,
    # then the one of g1 is used
    for attr_name_other in set(g2.vs.attributes()).intersection(set(g1.vs.attributes())).difference([v_attr_name]):
        attr_other = dict()
        for v in g2.vs():
            attr_other[attribute_to_ind[v[v_attr_name]]] = v[attr_name_other]

        for v in g1.vs():
            attr_other[attribute_to_ind[v[v_attr_name]]] = v[attr_name_other]
        graph_merged.vs[attr_name_other] = [attr_other[i] for i in range(graph_merged.vcount())]

    return graph_merged


def timed_main_run(main_function, logger_obj=logger):
    """Times the run of a function.

    :param main_function: The function to be run.
    :param logger_obj: The logger that will produce the result.
    :return:
    """
    start_time = time.time()
    logger_obj.info('Starting...')
    try:
        main_function()
    except Exception as e:
        logger_obj.error(type(e))
        logger_obj.error(traceback.format_exc())
    finally:
        logger_obj.info(f"Total time: {time.time() - start_time}")
