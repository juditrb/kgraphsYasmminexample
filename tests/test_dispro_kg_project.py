import os
import sys
sys.path.append("../")

import unittest
import pandas as pd

from biocypher import BioCypher
from dispro_kg_project import __version__
from dispro_kg_project.dispro_adapter import DisProAdapter

class TestLoad(unittest.TestCase):
    
    def setUp(self):
        all_data = pd.read_csv('../data/preprocessed_dispro.tsv', sep='\t')
        test_data = all_data.iloc[:100, :]
        test_data.to_csv("../data/testset.tsv", index=None, sep='\t')

        config_file = "../config/biocypher_config.yaml"
        data_file = "../data/testset.tsv"
        bc = BioCypher(
            biocypher_config_path = config_file,
        )
        self.adapter = DisProAdapter(
            data_file_path = data_file,
            biocypher_config_path = config_file,
            biocypher_driver = bc
        )
        
        ns = []
        for n in self.adapter.get_nodes():
            ns.append(n)
        self.nodes = ns
        
        es = []
        for e in self.adapter.get_edges():
            es.append(e)
        self.edges = es
    
    def tearDown(self):
        os.remove("../data/testset.tsv")    
        
    def test_version(self):
        self.assertTrue( __version__ == '0.1.0' )
       
    def test_quantity_nodes(self):
        self.assertTrue( len( self.nodes ) == 226 )
       
    def test_quantity_edges(self):
        self.assertTrue( len( self.edges ) == 313 )
        
    def test_node_existence(self):
        exnode = ('uniprot:P27695', 'protein', {'protein_name': 'DNA-(apurinic or apyrimidinic site) lyase'})
        self.assertTrue( exnode in self.nodes )
        
    def test_edge_existence(self):
        exedge = ('dispro:DP00007r002', 'dispro:DP00007', 'region_to_disorder', {'source': 'DisPro', 'version': 'v0.1', 'licence': 'MIT'})
        self.assertTrue( exedge in self.edges )
    
    def test_node_types(self):
        ntypes = {'disordered_region', 'organism', 'publication', 'evidence_type', 'disordered_peptide', 'biological_process', 'protein', 'structural_state', 'molecular_function', 'structural_transition', 'disorder_function'}
        ts = set( list(map( lambda x: x[1], self.nodes)) )
        
        self.assertTrue( len(ts) == 11 )
        self.assertTrue( len( ts.intersection(ntypes) ) == len(ntypes) )
    
    def test_edge_types(self):
        etypes = {'protein_to_region_disorder', 'region_to_disorder', 'protein_to_organism', 'region_to_evidence'}
        ts = set( list(map( lambda x: x[2], self.edges)) )
        
        self.assertTrue( len(ts) == 4 )
        self.assertTrue( len( ts.intersection(etypes) ) == len(etypes) )

if __name__ == '__main__':
    unittest.main()    
    
