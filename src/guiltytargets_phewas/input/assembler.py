# -*- coding: utf-8 -*-

"""Generates an adjacency file for the PPI network from the STRING database."""

import logging
from typing import Dict

import bio2bel_phewascatalog
import mygene
from opentargets import OpenTargetsClient
import pandas as pd
import psycopg2

from guiltytargets_phewas.constants import string_database, string_host, string_password, string_user
from guiltytargets_phewas.utils import get_converter_to_entrez, timed_main_run

logger = logging.getLogger(__name__)

__all__ = [
    'StringAssembler',
    'generate_phewas_file'
]


class StringAssembler:

    def __init__(self):
        self.ppi_url = f'https://stringdb-static.org/download/' \
            f'protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz'
        self.ppi_file = r'C:\Users\Mauricio\Thesis\data\raw\STRING\9606.protein.links.full.v11.0.txt.gz'

        self.out_cols = ['protein1', 'protein2', 'combined_score']

        # connection data to STRING database dump
        self.host = string_host
        self.database = string_database
        self.user = string_user
        self.password = string_password

    def create_adj_file(self, out_path, anno_type: str = 'entrezgene') -> None:
        """Create a csv file containing the protein interactions and the score.

        Requires the download from the protein network data for the organism Homo Sapiens from
        https://string-db.org/cgi/download.pl?sessionId=xoIrOw0ShV9o&species_text=Homo+sapiens.

        :param out_path: Path to the adjacency file (`string.edgelist`).
        :param anno_type: `entrezgene` for Entrez Id or `symbol` for Gene symbol.
        :return:
        """
        assert anno_type in ['symbol', 'entrezgene'], "STRING files can be created as either symbol or entrez ids"
        cols = ['protein1', 'protein2', 'neighborhood', 'neighborhood_transferred', 'fusion',
                'cooccurence', 'homology', 'coexpression', 'coexpression_transferred',
                'experiments', 'experiments_transferred', 'database', 'database_transferred',
                'textmining', 'textmining_transferred', 'combined_score']
        df = pd.read_csv(self.ppi_file, sep=' ', header=0, names=cols)

        ensp_to_symbol = self._get_ensembl_id_to_symbol_converter()

        prot1 = [ensp_to_symbol[x] for x in df['protein1']]
        prot2 = [ensp_to_symbol[x] for x in df['protein2']]

        if anno_type == 'entrezgene':
            unique_prots = set(prot1)
            unique_prots.union(prot2)
            symbol_to_entrez = get_converter_to_entrez(list(unique_prots))
            # not found symbols are not converted
            not_found = unique_prots.difference(symbol_to_entrez.keys())
            symbol_to_entrez.update({prot: prot for prot in not_found})
            prot1 = [symbol_to_entrez[x] for x in prot1]
            prot2 = [symbol_to_entrez[x] for x in prot2]

        df.loc[:, 'prot_symbol1'] = pd.Series(prot1)
        df.loc[:, 'prot_symbol2'] = pd.Series(prot2)

        df[self.out_cols].to_csv(out_path, sep='\t', header=False, index=False)

    def _get_ensembl_id_to_symbol_converter(self) -> Dict[str, str]:
        """Get a converter from Ensembl Id (ENSPxxx) to protein symbol.

        Requires previous installation of PostGreSQL and to import the database schema ``items''
        from STRING for organism Homo Sapiens (download at
        https://string-db.org/cgi/download.pl?sessionId=xoIrOw0ShV9o&species_text=Homo+sapiens).

        :return: a dictionary with keys as Ensembl identifiers and values as protein symbols.
        """
        schema = 'items'
        sql_command = " SELECT protein_external_id, protein_id, preferred_name"
        sql_command += "  FROM {}.{}".format(str(schema), str('proteins'))
        sql_command += ' WHERE species_id = {}'.format(str(9606))
        sql_command += ';'

        conn = psycopg2.connect(
            host=self.host,
            database=self.database,
            user=self.user,
            password=self.password
        )
        data = pd.read_sql(sql_command, conn)
        conn.close()

        return dict(data[['protein_external_id', 'preferred_name']].values)


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

    # TODO use the converters.get_converter_to_entrez
    mg = mygene.MyGeneInfo()
    id_mappings = mg.getgenes(ensembl_list, fields=anno_type)

    with open(outpath, 'w+') as outfile:
        for mapping in id_mappings:
            if anno_type in mapping.keys():
                outfile.write(mapping[anno_type])
                outfile.write('\n')


def main():
    assembler = StringAssembler()
    assembler.create_adj_file('string_entrez1.edgelist', 'entrezgene')
    # input.create_adj_file('string_symbol.edgelist', 'symbol')


if __name__ == '__main__':
    timed_main_run(main, logger)
