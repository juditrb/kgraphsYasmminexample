
import pandas as pd
from bioregistry import normalize_curie

def run():
    df = pd.read_csv('../data/dispro_release.tsv', sep='\t')

    annot = { "BIOLOGICAL_PROCESS": "biological_process", "CELLULAR_COMPONENT": "cellular_component", "MOLECULAR_FUNCTION": "molecular_function", "STRUCTURAL_STATE": "structural_state", "STRUCTURAL_TRANSITION": "structural_transition", "DISORDER_FUNCTION": "disorder_function" }
    
    namespaces = { "acc": "uniprot", "ncbi_taxon_id": "ncbitaxon", "disprot_id": "dispro", "region_id": "dispro" }
    norm = ["ec", "reference"]
    
    aux = {}
    for a in annot:
        for fs in ['id', 'name']:
            aux[ f"{annot[a]}-{fs}" ] = []
    for a in namespaces:
        aux[ a ] = []
    for a in norm:
        aux[ a ] = []
        
    types_annot = df['term_namespace'].tolist()

    for i in df.index:
        id = df.loc[i, 'term']
        name = df.loc[i, 'term_name']
        ann = annot[ types_annot[i].upper().replace(' ','_') ]
        
        for a in annot:
            for fs in ['id', 'name']:
                if( ann == annot[a] ):
                    v = eval(fs)
                    if( fs=='id' ):
                        v = normalize_curie( v )
                    aux[ f"{ann}-{fs}" ].append( v ) 
                else:
                    aux[ f"{annot[a]}-{fs}" ].append( "" ) 
        
        for a in namespaces:
            v = f"{ namespaces[a] }:{ df.loc[i, a] }"
            v = normalize_curie(v)
            if( v == None ):
                v = f"{ namespaces[a] }:{ df.loc[i, a] }"
            aux[ a ].append( v )
            
        for a in norm:
            v = df.loc[i, a]
            v = normalize_curie(v)
            aux[ a ].append( v )
        
    filt = df[ ['name', 'organism', 'start', 'end', 'ec_name', 'region_sequence'] ]
    filt['name'] = filt.name.apply(lambda x: x.replace("'", ''))
    for k in aux:
        filt[k] = aux[k]
    filt.to_csv( '../data/preprocessed_dispro.tsv', index=None, sep='\t' )
    
run()
