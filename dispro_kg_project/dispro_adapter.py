#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - Disordered Proteins adapter prototype
"""

import csv
import yaml
import neo4j_utils as nu
import pandas as pd

from enum import Enum
from typing import Optional
from bioregistry import normalize_curie
from itertools import chain

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class DisProNodeType(Enum):
    """
    DisPro nodes.
    """

    PROTEIN = "protein"
    ORGANISM = "organism"
    EVIDENCE_TYPE = "evidence_type"
    PUBLICATION = "publication"
    DISORDER_PEPTIDE = "disordered_peptide"
    DISORDER_REGION = "disordered_region"
    BIOLOGICAL_PROCESS = "biological_process"
    CELLULAR_COMPONENT = "cellular_component"
    MOLECULAR_FUNCTION = "molecular_function"
    STRUCTURAL_STATE = "structural_state"
    STRUCTURAL_TRANSITION = "structural_transition"
    DISORDER_FUNCTION = "disorder_function"

class DisProProteinNodeField(Enum):
    """
    Fields available for DisPro proteins.
    """

    ID = "acc"
    PROTEIN_NAME = "name"

class DisProOrganismNodeField(Enum):
    """
    Fields available for NCBI Taxon organism.
    """

    ID = "ncbi_taxon_id"
    ORGANISM_NAME = "organism"
    
class DisProEvidenceTypeNodeField(Enum):
    """
    Fields available for Evidences.
    """

    ID = "ec"
    EVIDENCE_NAME = "ec_name"
    
class DisProPublicationNodeField(Enum):
    """
    Fields available for Evidences.
    """

    ID = "reference"

class DisProDisorderPeptideNodeField(Enum):
    """
    Fields available for Disordered Peptides.
    """

    ID = "disprot_id"

class DisProDisorderRegionNodeField(Enum):
    """
    Fields available for Disordered Regions.
    """

    ID = "region_id"
    DISORDER_REGION_START = "start"
    DISORDER_REGION_END = "end"
    DISORDER_REGION_SEQUENCE = "region_sequence"
    
class DisProGOBPNodeField(Enum):
    """
    Fields available for GO BP.
    """

    ID = "biological_process-id"
    ANNOT_BP_NAME = "biological_process-name"
    
class DisProGOCCNodeField(Enum):
    """
    Fields available for GO CC.
    """

    ID = "cellular_component-id"
    ANNOT_CC_NAME = "cellular_component-name"
    
class DisProGOMFNodeField(Enum):
    """
    Fields available for GO MF.
    """

    ID = "molecular_function-id"
    ANNOT_mf_NAME = "molecular_function-name"
    
class DisProIDPOSSNodeField(Enum):
    """
    Fields available for IDPO Structural State.
    """

    ID = "structural_state-id"
    ANNOT_ss_NAME = "structural_state-name"
    
class DisProIDPOSTNodeField(Enum):
    """
    Fields available for IDPO Structural Transition.
    """

    ID = "structural_transition-id"
    ANNOT_ST_NAME = "structural_transition-name"
    
class DisProIDPODFNodeField(Enum):
    """
    Fields available for IDPO Disorder Function.
    """

    ID = "disorder_function-id"
    ANNOT_DF_NAME = "disorder_function-name"

class DisProEdgeType(Enum):
    """
    DisPro edges.
    """

    PROTEIN_TO_ORGANISM = "protein_to_organism"
    REGION_TO_DISORDER = "region_to_disorder"
    REGION_TO_EVIDENCE = "region_to_evidence"
    PROTEIN_TO_REGION_DISORDER = "protein_to_region_disorder"

class DisProProteinToOrganismEdgeField(Enum):
    """
    Fields available for DisPro protein to organism edges.
    """
    
    SOURCE = "acc"
    TARGET = "ncbi_taxon_id"

class DisProRegionToDisorderEdgeField(Enum):
    """
    Fields available for DisPro region to disorder edges.
    """
    
    SOURCE = "region_id"
    TARGET = "disprot_id"

class DisProRegionToEvidenceEdgeField(Enum):
    """
    Fields available for DisPro region to evidence edges.
    """
    
    SOURCE = "region_id"
    TARGET = "ec"
    
class DisProProteinToRegionDisorderEdgeField(Enum):
    """
    Fields available for DisPro protein to region edges.
    """

    SOURCE = "acc"
    TARGET = "region_id"
    REFERENCE_PUBLICATION = "reference"
    FUNCTIONAL_ANNOTATION = "annotation_term"

class DisProAdapter:
    def __init__(
        self,
        data_file_path: str,
        biocypher_config_path: str,
        biocypher_driver = None,
        #schema_file_path: str,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        # 'data/preprocessed_dispro.tsv'
        # 'config/schema_config.yml'
        self.data = pd.read_csv( data_file_path, sep='\t')
        #self.schema = yaml.safe_load( open( schema_file_path,'r') )
        
        self._set_up_types_and_fields(
            node_types, node_fields, edge_types, edge_fields
        )
        
        self.annotation_node_fields = [ DisProGOBPNodeField, DisProGOCCNodeField, DisProGOMFNodeField, DisProIDPOSSNodeField, DisProIDPOSTNodeField, DisProIDPODFNodeField ]
        
        self.biocypher_driver = biocypher_driver
        
        if( not self.biocypher_driver.base_config['offline'] ):
            self.dbconfig = yaml.safe_load( open( biocypher_config_path,'r') )['neo4j']
            
            # read driver
            self.driver = nu.Driver(
                db_name = self.dbconfig.database_name,
                db_uri = self.dbconfig.uri,
                db_user = self.dbconfig.user,
                db_passwd = self.dbconfig.password,
                multi_db = False,
                max_connection_lifetime = 7200,
            )
            if( not self.db_exists( self.dbconfig.database_name ) ):
                self.create_db( self.dbconfig.database_name )
            
            else:
                if( self.dbconfig.wipe ):
                    self.driver.wipe_db()
            
            self.select_db( self.dbconfig.database_name )
            
        self.data_source = "DisPro"
        self.data_version = "v0.1"
        self.data_licence = "MIT"

    def _parse_info_node(self, row, fields):
        _id = row[ fields['ID'].value ]
        _properties = {}
        for f in fields:
            if( (not f.name in ['ID']) and (f.value in self.node_fields) ):
                key = f.name.lower()
                value = row[ f.value ]
                if( str(value) != 'nan' ):
                    if( str(value).find("'") != -1 ):
                        value = value.replace("'",'')
                    _properties[ key ] = value
        
        return _id, _properties
    
    def get_nodes(self):
        """
        Get nodes from CSV and yield them to the batch writer.

        Args:
            label: input label of nodes to be read

        Returns:
            generator of tuples representing nodes
        """

        loc_dict = {
            'protein': DisProProteinNodeField,
            'organism': DisProOrganismNodeField,
            'evidence_type': DisProEvidenceTypeNodeField,
            'publication': DisProPublicationNodeField,
            'disordered_peptide': DisProDisorderPeptideNodeField,
            'disordered_region': DisProDisorderRegionNodeField,
            'biological_process': DisProGOBPNodeField,
            'cellular_component': DisProGOCCNodeField,
            'molecular_function': DisProGOMFNodeField,
            'structural_state': DisProIDPOSSNodeField,
            'structural_transition': DisProIDPOSTNodeField,
            'disorder_function': DisProIDPODFNodeField
        }

        for ntype in self.node_types:
            nodes = set()
            fields = loc_dict[ ntype ]
            for index, row in self.data.iterrows():
                _id, _props = self._parse_info_node(row, fields)
                if( not _id in nodes ):
                    nodes.add(_id)
                    if( str(_id) != 'nan' ):
                        yield _id, ntype, _props

    def _get_annotation_nodeId(self, row):
        _id = None
        for f in self.annotation_node_fields:
            _id = row[ f.ID.value ]
            if( (_id != '') and (str(_id) != 'nan') and (_id != None) ):
                return _id
        
        return _id        
    
    def _parse_info_edge(self, row, fields):
        _source = row[ fields['SOURCE'].value ]
        _target = row[ fields['TARGET'].value ]
        
        _properties = {}
        for f in fields:
            if( (not f.name in ['SOURCE', 'TARGET']) and (f.value in self.edge_fields) ):
                key = f.name
                if( key == 'FUNCTIONAL_ANNOTATION' ):
                    value = self._get_annotation_nodeId(row)
                else:
                    value = row[ f.value ]
                _properties[ key.lower() ] = value
        
        _properties["source"] = self.data_source
        _properties["version"] = self.data_version
        _properties["licence"] = self.data_licence
        
        return _source, _target, _properties

    def get_edges(self):
        """
        Get edges from CSV and yield them to the batch writer.

        Args:
            label: input label of edges to be read

        Returns:
            generator of tuples representing edges
        """

        loc_dict = {
            'protein_to_organism': DisProProteinToOrganismEdgeField,
            'region_to_disorder': DisProRegionToDisorderEdgeField,
            'region_to_evidence': DisProRegionToEvidenceEdgeField,
            'protein_to_region_disorder': DisProProteinToRegionDisorderEdgeField
        }

        for etype in self.edge_types:
            edges = set()
            _label = etype
            fields = loc_dict[ _label ]
            for index, row in self.data.iterrows():
                _source, _target, _props = self._parse_info_edge(row, fields)
                _id = f"{_source}_{_target}"
                if( not _id in edges ):
                    edges.add(_id)
                    yield _source, _target, _label, _props
    
    def _set_up_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = [field.value for field in node_types]
        else:
            self.node_types = [field.value for field in DisProNodeType]

        if edge_types:
            self.edge_types = [field.value for field in edge_types]
        else:
            self.edge_types = [field.value for field in DisProEdgeType]

        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    DisProProteinNodeField,
                    DisProOrganismNodeField,
                    DisProEvidenceTypeNodeField,
                    DisProPublicationNodeField,
                    DisProDisorderPeptideNodeField,
                    DisProDisorderRegionNodeField,
                    DisProGOBPNodeField,
                    DisProGOCCNodeField,
                    DisProGOMFNodeField,
                    DisProIDPOSSNodeField,
                    DisProIDPOSTNodeField,
                    DisProIDPODFNodeField
                )
            ]
        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [
                field.value
                for field in chain(
                    DisProProteinToOrganismEdgeField,
                    DisProRegionToDisorderEdgeField,
                    DisProRegionToEvidenceEdgeField,
                    DisProProteinToRegionDisorderEdgeField
                )
            ]

# multi-line fields: only due to line 832 in cellModels_all.csv?
