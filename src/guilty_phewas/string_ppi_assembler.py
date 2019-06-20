# -*- coding: utf-8 -*-

"""Generates an adjacency file for the PPI network from the STRING database."""

from typing import Dict
import pandas as pd
import psycopg2


class StringAssembler:

    def __init__(self):
        self.adj_file = 'string.edgelist'
        self.ppi_url = f'https://stringdb-static.org/download/' \
            f'protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz'
        self.ppi_file = 'C:\\Users\\Mauricio\\Thesis\\9606.protein.links.full.v11.0.txt.gz'

        self.out_cols = ['prot_symbol1', 'prot_symbol2', 'combined_score']

        self.conn = psycopg2.connect(
            host="localhost",
            database="stringitems",
            user="postgres",
            password="123456"
        )

        self._dir = ""

    def create_adj_file(self):
        """
        Creates a csv file containing the protein interactions and the score.
        Requires the download from the protein network data for the organism Homo Sapiens from
        https://string-db.org/cgi/download.pl?sessionId=xoIrOw0ShV9o&species_text=Homo+sapiens

        :param url:
        :param out_path:
        :return:
        """
        # s = requests.get(url).content

        cols = ['protein1', 'protein2', 'neighborhood', 'neighborhood_transferred', 'fusion',
                'cooccurence', 'homology', 'coexpression', 'coexpression_transferred',
                'experiments', 'experiments_transferred', 'database', 'database_transferred',
                'textmining', 'textmining_transferred', 'combined_score']
        df = pd.read_csv(self.ppi_file, sep=' ', header=0, names=cols)
        ensp_to_symbol = self._get_ensembl_id_to_symbol_converter()
        df.loc[:, 'prot_symbol1'] = pd.Series([ensp_to_symbol[x] for x in df['protein1']])
        df.loc[:, 'prot_symbol2'] = pd.Series([ensp_to_symbol[x] for x in df['protein2']])
        df[self.out_cols].to_csv(self.adj_file, sep='\t', header=False, index=False)

    def _get_ensembl_id_to_symbol_converter(self) -> Dict:
        """
        Get a converter from Ensembl Id (ENSPxxx) to protein symbol. Requires previous
        installation of PostGreSQL and to import the database schema ``items''
        from STRING for organism Homo Sapiens (download at
        https://string-db.org/cgi/download.pl?sessionId=xoIrOw0ShV9o&species_text=Homo+sapiens).

        :return: a dictionary with keys as Ensembl identifiers and values as protein symbols.
        """
        schema = 'items'
        sql_command = " SELECT protein_external_id, protein_id, preferred_name"
        sql_command += "  FROM {}.{}".format(str(schema), str('proteins'))
        sql_command += ' WHERE species_id = {}'.format(str(9606))
        sql_command += ';'

        data3 = pd.read_sql(sql_command, self.conn)

        return {x: y for x, y in data3[['protein_external_id', 'preferred_name']].values}


if __name__ == '__main__':
    assembler = StringAssembler()
    assembler.create_adj_file()
