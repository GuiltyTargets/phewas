# -*- coding: utf-8 -*-

"""Utility functions to run the algorithm"""

from typing import Dict, Mapping, Optional

from networkx import DiGraph
from pybel.dsl import gene, protein, rna


def add_disease_attribute(graph: DiGraph, att_mappings: Dict):
    """Add the phenotypes to the Base Entities as attributes."""
    for node in graph.nodes:
        if ((isinstance(node, protein) or
             isinstance(node, rna) or
             isinstance(node, gene)) and
                node.name in att_mappings):
            graph.nodes[node]['phenotypes'] = [phtype for _, phtype in att_mappings[node.name]]


def write_adj_file_attribute(graph, filepath: str, att_mappings: Dict):
    """Write an adjacency file from the attribute graph."""
    with open(filepath, 'w') as f:
        for node in graph.nodes:
            if 'phenotypes' in graph.nodes[node]:  # "There are diseases in the node":
                print(f"{node} {' '.join(str(att_mappings[phe]) for phe in graph.nodes[node]['phenotypes'])}", file=f)
