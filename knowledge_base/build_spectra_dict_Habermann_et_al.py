# Filippo Gastaldello - 11/07/2026

# Build gene set dictionary for Spectra using cytopus

import cytopus as cp
import networkx as nx
import matplotlib.pyplot as plt
import json

# Save path for spectra dictionary and conversion dictionary
spectra_dict_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Habermann_et_al/cytopus_gene_sets/knowledgebase_gene_sets_manual_annotation.json"
conversion_dict_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Habermann_et_al/cell_type_conversion_dictionaries/Habermann_et_al_conversion_dict_mixed_annotation_to_cytopus.json"

# Load KnowledgeBase
G = cp.KnowledgeBase()

# Build conversion dictionary from the dataset's cell types to the knowledgebase's cell types
conversion_dict = {
    "MUC5AC+ High" : "lung-epi",                  
    "Basal" : "lung-epi",    
    "Proliferating Epithelial Cells" : "lung-epi",
    "Differentiating Ciliated" : "lung-epi",      
    "Ciliated" : "lung-epi",                      
    "SCGB3A2+ SCGB1A1+" : "lung-epi",             
    "SCGB3A2+" : "lung-epi",                      
    "MUC5B+" : "lung-epi",                        
    "KRT5-/KRT17+" : "lung-epi",                  
    "Transitional AT2" : "lung-epi",                  
    "AT2" : "lung-epi",                           
    "AT1" : "lung-epi",                           
    "Lymphatic Endothelial Cells" : "endo-lymphatic",   
    "Endothelial Cells" : "endo",             
    "Proliferating Macrophages" : "Mac",     
    "Proliferating T Cells" : "T",         
    "NK Cells" : "NK",  
    "pDCs" : "p-DC",                          
    "cDCs" : "cDC",                          
    "Mast Cells" : "mast",                    
    "Plasma Cells" : "plasma",                  
    "B Cells" : "B",                       
    "Monocytes" : "mono",                     
    "T Cells" : "T",                       
    "Macrophages" : "Mac",                   
    "PLIN2+ Fibroblasts" : "fibro",            
    "Fibroblasts" : "fibro",                   
    "HAS1 High Fibroblasts" : "fibro",         
    "Myofibroblasts" : "fibro",                
    "Mesothelial Cells" : "lung-epi",             
    "Smooth Muscle Cells" : "lung-smooth-muscle"
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
