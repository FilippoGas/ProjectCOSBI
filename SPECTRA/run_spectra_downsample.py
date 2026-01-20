import numpy as np
import json 
import scanpy as sc
from collections import OrderedDict
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import csv
import pickle
import h5py
import Spectra as spc
import sys
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets


def main():

    # SCRIPT VARIABLES
    downsampling        = sys.argv[1]
    annotation_type     = sys.argv[2]
    annotation_column   = sys.argv[3]
    pbs_job_id          = sys.argv[4]
    # INPUT PATHS
    scKnowledgeBase_dict_path   = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cytopus_gene_sets/knowledgebase_gene_sets_"+annotation_type+"_annotation.json"
    sc_dataset_path             = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/AnnData/Natri_et_al_downsampled_"+downsampling+".h5ad"

    # OUTPUT PATHS
    model_pickle_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_"+downsampling+"/"+annotation_type+"_annotation/Natri_et_al_downsampled_"+downsampling+"_"+annotation_type+"_anno_spectra_model_"+pbs_job_id+".pickle"
    anndata_pickle_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_"+downsampling+"/"+annotation_type+"_annotation/Natri_et_al_downsampled_"+downsampling+"_"+annotation_type+"_anno_AnnData_"+pbs_job_id+".pickle"
    gene_scores_csv_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_"+downsampling+"/"+annotation_type+"_annotation/Natri_et_al_downsampled_"+downsampling+"_"+annotation_type+"_anno_gene_scores_"+pbs_job_id+".csv"
    cell_scores_csv_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_"+downsampling+"/"+annotation_type+"_annotation/Natri_et_al_downsampled_"+downsampling+"_"+annotation_type+"_anno_cell_scores_"+pbs_job_id+".csv"
    factor_markers_csv_path = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_"+downsampling+"/"+annotation_type+"_annotation/Natri_et_al_downsampled_"+downsampling+"_"+annotation_type+"_anno_factor_markers_"+pbs_job_id+".csv"

    # LOAD DATA
    # Load gene sets dictionary from cytopus
    with open(scKnowledgeBase_dict_path, 'r') as f:
        scKnowledgeBase_dict = json.load(f)

    # Before reading the h5ad file I need to remove the attribute "active.ident" added by SeuratDisk during the save of the seurat object into a h5ad file
    # (might not be needed but keep this for the moment)
    with h5py.File(sc_dataset_path, "r") as f:
        # Check if 'active.ident' exists as file attribute and delete it
        if "active.ident" in f.attrs:
            del f.attrs["active.ident"]
    
    # Read scRNAseq dataset as AnnData
    adata = sc.read_h5ad(sc_dataset_path)
    
    # FIT SPECTRA MODEL
    print("\nStart model training....\n") 
    trained_model = spc.est_spectra(adata=adata,
                                    gene_set_dictionary=scKnowledgeBase_dict,
                                    L=None,
                                    use_highly_variable=False,
                                    cell_type_key=annotation_column,
                                    use_weights=True,
                                    lam=0.1,
                                    delta=0.001,
                                    kappa=None,
                                    rho=0.001,
                                    use_cell_types=True,
                                    n_top_vals=50,
                                    label_factors=True,
                                    overlap_threshold=0.2,
                                    clean_gs=True,
                                    min_gs_num=3,
                                    num_epochs=10000) 
    print("\nDone model training.\n")

    # SAVE RAW RESULTS
    # Save trained model
    with open(model_pickle_path, "wb") as f:
        pickle.dump(trained_model, f)
    # Save adata
    with open(anndata_pickle_path, "wb") as f:
        pickle.dump(adata, f)

    # Save result matrices
    gene_scores     =   pd.DataFrame(adata.uns['SPECTRA_factors'], 
                                 index = adata.uns['SPECTRA_overlap'].index, 
                                 columns = list(adata.var.axes[0]))             # factors x genes matrix that tells you how important each gene is to the resulting factors 
    markers     =   pd.DataFrame(adata.uns['SPECTRA_markers'].T, 
                                 columns = adata.uns['SPECTRA_overlap'].index)  # factors x n_top_vals list of n_top_vals top markers per factor
    cell_scores =   pd.DataFrame(adata.obsm['SPECTRA_cell_scores'], 
                                 index = adata.obs_names, 
                                 columns = adata.uns['SPECTRA_overlap'].index)  # cells x factors matrix of cell scores
    gene_scores.to_csv(gene_scores_csv_path)
    markers.to_csv(factor_markers_csv_path)
    cell_scores.to_csv(cell_scores_csv_path)

if __name__ == "__main__":
    main()
