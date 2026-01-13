# Filippo Gastaldello - 19/11/2025

# Build gene set dictionary for Spectra using cytopus

import cytopus as cp
import networkx as nx
import matplotlib.pyplot as plt
import json

# Save path for spectra dictionary and conversion dictionary
spectra_dict_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cytopus_gene_sets/knowledgebase_gene_sets_manual_annotation.json"
conversion_dict_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cell_type_conversion_dictionaries/Natri_et_al_conversion_dict_manual_annotation_to_cytopus.json"

# Load KnowledgeBase
G = cp.KnowledgeBase()

# Build conversion dictionary from the dataset's cell types to the knowledgebase's cell types
conversion_dict = {
    'AT2':'lung-epi',
    'Proliferating - Epi':'lung-epi',
    'Ciliated':'lung-epi',
    'AT1':'lung-epi',
    'Secretory - SCGB3A2+':'lung-epi',
    'Basal':'lung-epi',
    'Transitional AT2':'lung-epi',
    'Secretory - SCGB1A1+/MUC5B+':'lung-epi',
    'Secretory - SCGB1A1+/SCGB3A2+':'lung-epi',
    'Differentiating ciliated':'lung-epi',
    'KRT5-/KRT17+':'lung-epi',
    'PNEC':'lung-epi',
    'PLIN2+ FB':'fibro',
    'MyoFB':'fibro',
    'Adventitial FB':'fibro',
    'Alveolar FB':'fibro',
    'SMC':'lung-smooth-muscle',
    'Pericyte':'smooth-muscle',
    'Mesothelial':'epi',
    'MyoFB - Activated':'fibro',
    'Venule':'lung-endo-venous',
    'Arteriole':'endo-arterial',
    'Lymphatic':'endo-lymphatic',
    'Systemic venous':'endo-systemic-venous',
    'aCap':'capillary',
    'gCap':'capillary',
    'Inflamed':'leukocyte',
    'Monocyte':'mono',
    'Monocyte-derived macrophage':'mono',
    'NK':'NK',
    'Inflammatory monocyte':'mono',
    'moDC':'mo-DC',
    'B cells':'B',
    'cDC1':'cDC1',
    'cDC2':'cDC2',
    'Alveolar macrophage':'Mac',
    'Proliferating - Imm':'leukocyte',
    'Plasma':'plasma',
    'Interstitial macrophage':'Mac',
    'CD4':'CD4-T',
    'Mast':'mast',
    'CD8/NKT':'T',
    'pDC':'p-DC'
}
# Save conversion dict
with open(conversion_dict_path, 'w') as f:
    json.dump(conversion_dict, f)

# Get cell types present in the dataset
celltype_of_interest = list(dict.fromkeys(conversion_dict.values()))
global_celltype = ['all-cells']

# Get spectra dictionary
G.get_celltype_processes(celltype_of_interest, global_celltypes=global_celltype)
processes_dict = G.celltype_process_dict

# Convert lables back to the original version and save spectra dict
spectra_dict = {}
for key, value in conversion_dict.items():
    spectra_dict[key] = processes_dict[value]
spectra_dict["global"] = processes_dict["global"]
with open(spectra_dict_path, 'w') as f:
    json.dump(spectra_dict, f)
