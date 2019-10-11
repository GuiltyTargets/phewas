# -*- coding: utf-8 -*-

""""""

import logging
import mygene
from typing import Dict, List, Set, Union
from urllib import parse, request

logger = logging.getLogger(__name__)


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