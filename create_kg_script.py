import sys
sys.path.append("")

import os
import io
import pstats
import cProfile

from dispro_kg_project.dispro_adapter import DisProAdapter
from biocypher import BioCypher

def main():
    """
    Connect BioCypher to DepMap adapter to import data into Neo4j.

    Optionally, run with profiling.
    """
    
    PROFILE = True
    
    if PROFILE:
        profile = cProfile.Profile()
        profile.enable()
        
    ###############
    # ACTUAL CODE #
    ###############
    
    cwd = os.getcwd()
    config_file = os.path.join( cwd, "config/biocypher_config.yaml" )
    schema_file = os.path.join( cwd, "config/schema_config.yaml" )
    data_file = os.path.join( cwd, "data/preprocessed_dispro.tsv" )
    
    # start biocypher
    bc = BioCypher(
        biocypher_config_path = config_file,
    )

    # check schema
    #bc.show_ontology_structure() # Eco is too big

    # create adapter
    dispro = DisProAdapter(
        data_file_path = data_file,
        biocypher_config_path = config_file,
        biocypher_driver = bc
    )

    # write nodes and edges to csv
    bc.write_nodes( dispro.get_nodes() )
    bc.write_edges( dispro.get_edges() )

    # convenience and stats
    bc.write_schema_info(as_node=True)
    bc.write_import_call()
    
    bc.log_missing_input_labels()
    bc.log_duplicates()
    # bc.summary() It shows the ontology structure too

    ######################
    # END OF ACTUAL CODE #
    ######################

    if PROFILE:
        profile.disable()

        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(profile, stream=s).sort_stats(sortby)
        ps.print_stats()

        ps.dump_stats("adapter.prof")
        # look at stats using snakeviz


if __name__ == "__main__":
    main()
    
