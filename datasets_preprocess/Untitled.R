
library(tidyverse)
library(Seurat)

data_dir = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Habermann_et_al/Habermann_et_al.rds"
out_dir <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Habermann_et_al/"


seurat_obj <- readRDS("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Habermann_et_al/Habermann_et_al.rds")

tokeep <- seurat_obj_metadata %>%
        rownames_to_column(var = "UMI") %>%
        dplyr::filter(Diagnosis %in% c("Control","IPF") & percent.mt < 20) %>%
        dplyr::pull(UMI)


print(seurat_obj)
seurat_obj <- seurat_obj[,colnames(seurat_obj) %in% tokeep]
print(seurat_obj)


#normalize and scale
#seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
#seurat_obj <- ScaleData(seurat_obj,vars.to.regress = c("percent.mt","nCount_RNA"))

print(names(seurat_obj@meta.data))

# sctrasform of the subsetted object
seurat_obj <- SCTransform(seurat_obj, batch_var = "orig.ident",method = 'glmGamPoi',vars.to.regress = c("percent.mt","nCount_RNA"))


#save metadata
saveRDS(seurat_obj@meta.data, file=paste0(out_dir,"seurat_GSE135893_filtered_metadata.rds"))
#save filtered seurat object
saveRDS(seurat_obj, file=paste0(out_dir,"seurat_GSE135893_filtered.rds"))
