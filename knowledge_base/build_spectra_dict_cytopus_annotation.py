# Filippo Gastaldello - 19/11/2025

# Build gene set dictionary for Spectra using cytopus

import cytopus as cp
import networkx as nx
import matplotlib.pyplot as plt
import json
import pandas as pd

# Save path for spectra dictionary and conversion dictionary
spectra_dict_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cytopus_gene_sets/knowledgebase_gene_sets_cytopus_annotation.json"

# Load KnowledgeBase
G = cp.KnowledgeBase()

# Get cell types present in the dataset with cytopus annotation
celltype_of_interest = list(set(pd.read_csv("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cell_type_lists/Natri_et_al_annotated_cell_types_qc.csv")["cytopus"]))
global_celltype = ['all-cells']

# Get spectra dictionary
G.get_celltype_processes(celltype_of_interest, global_celltypes=global_celltype)
spectra_dict = G.celltype_process_dict

# Save spectra dict
with open(spectra_dict_path, 'w') as f:
    json.dump(spectra_dict, f)
