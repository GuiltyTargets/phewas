# -*- coding: utf-8 -*-

"""Parser methods for the input files."""

import logging
import os

from igraph import Graph
import pandas as pd
from pybel import BELGraph, union

from .converters import get_converter_to_disease_ontology, uniprot_id_to_entrez_id_converter

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)


def parse_disease_gene_graph(path: str) -> Graph:
    """ """
    logger.info("Loading Disease Gene Graph.")
    DISEASE_COLUMN = '# Disease ID'
    GENE_COLUMN = 'Gene ID'
    data = pd.read_csv(path, sep='\t', dtype=str)
    data[DISEASE_COLUMN] = [
        'UMLS_CUI:' + disease
        for disease
        in data[DISEASE_COLUMN]
    ]
    doid_conveter = get_converter_to_disease_ontology()
    data[DISEASE_COLUMN] = [
        doid_conveter[disease] if disease in doid_conveter else disease
        for disease
        in data[DISEASE_COLUMN]
    ]
    logger.debug(f"# of diseases doid'ed: {len({x for x in data[DISEASE_COLUMN] if x.startswith('DOID:')})}")
    logger.debug(f"# of not doid'ed:      {len({x for x in data[DISEASE_COLUMN] if x.startswith('UMLS_CUI:')})}")
    columns = [DISEASE_COLUMN, GENE_COLUMN]

    edgelist = [list(edge) for _, edge in data[columns].iterrows()]
    edgelist.extend([[edge[1], edge[0]] for _, edge in data[columns].iterrows()])
    return Graph.TupleList(edges=edgelist, directed=False)


def parse_disease_drug_graph(path: str) -> Graph:
    """ """
    logger.info("Loading Disease Drug Graph.")
    DISEASE_COLUMN = '# Disease(MESH)'
    DRUG_COLUMN = 'Chemical'
    data = pd.read_csv(path, sep='\t', dtype=str)
    doid_conveter = get_converter_to_disease_ontology()
    data[DISEASE_COLUMN] = [
        doid_conveter[disease] if disease in doid_conveter else disease
        for disease
        in data[DISEASE_COLUMN]
    ]
    columns = [DISEASE_COLUMN, DRUG_COLUMN]
    logger.debug(f"# of diseases doid'ed: {len({x for x in data[DISEASE_COLUMN] if x.startswith('DOID:')})}")
    logger.debug(f"# of not doid'ed:      {len({x for x in data[DISEASE_COLUMN] if x.startswith('MESH:')})}")

    edgelist = [list(edge) for _, edge in data[columns].iterrows()]
    edgelist.extend([[edge[1], edge[0]] for _, edge in data[columns].iterrows()])
    return Graph.TupleList(edges=edgelist, directed=False)


def parse_gene_drug_graph(path: str) -> Graph:
    """ """
    logger.info("Loading Drug Gene Graph.")
    DRUG_COLUMN = '#Drug'
    GENE_COLUMN = 'Gene'
    data = pd.read_csv(path, sep='\t', dtype=str)
    columns = [GENE_COLUMN, DRUG_COLUMN]

    uniprot_ids = set(data[GENE_COLUMN])
    converter = uniprot_id_to_entrez_id_converter(uniprot_ids)
    data[GENE_COLUMN] = [
        converter[x] if x in converter else x
        for x in data[GENE_COLUMN]
    ]

    edgelist = [list(edge) for _, edge in data[columns].iterrows()]
    edgelist.extend([[edge[1], edge[0]] for _, edge in data[columns].iterrows()])
    return Graph.TupleList(edges=edgelist, directed=False)


def bel_graph_loader(from_dir: str) -> BELGraph:
    """Obtains a combined BELGraph from all the BEL documents in one folder.

    :param from_dir: The folder with the BEL documents.
    :return: A corresponding BEL Graph.
    """
    logger.info("Loading BEL Graph.")
    files = [
        os.path.join(from_dir, file)
        for file
        in os.path.listdir(from_dir)
        if os.path.isfile(os.path.join(from_dir, file))
    ]
    bel_files = [file for file in files if file[-4:].lower() == '.bel']
    bel_graphs = [os.path.from_path(file) for file in bel_files]
    return union(bel_graphs)
