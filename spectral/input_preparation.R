# Filippo Gastaldello - 20/01/2026

# Prepare input for spektral - Takes SPECTRA output and divide it into one table
# per patient, at cell type resolution

library(tidyverse)
library(Seurat)         # To retrieve cell type annotations
library(jsonlite)

# Load SPECTRA output
cell_scores <- read_csv("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/downsampling_4000/cytopus_annotation/Natri_et_al_downsampled_4000_cytopus_anno_cell_scores_14232.hpc3-head1.unitn.it.csv")
colnames(cell_scores)[1] <- "cell_id"
# Load seurat object for cell type annotations and cell-patient match
data <- readRDS("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/Natri_et_al.rds")

# Outpout folder
output_folder <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/spektral/inputs/"
dir.create(paste0(output_folder, "patients_cell_scores"), recursive = TRUE)

# Add annotation from cytopus
cytopus_dict <- read_json("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cell_type_conversion_dictionaries/Natri_et_al_conversion_dict_manual_annotation_to_cytopus.json")
cytopus_dict <- data.frame("manual_annotation_1"=names(cytopus_dict),
                           "cytopus"=unname(unlist(cytopus_dict))
                           )
data@meta.data <- data@meta.data %>% 
                        rownames_to_column("barcode") %>%
                        left_join(cytopus_dict, by = "manual_annotation_1") %>%
                        column_to_rownames("barcode")
# Extract metadata and only keep control and IPF cells
metadata <- data@meta.data %>%
                select(Sample_Name, Diagnosis, Sample_Type, Tobacco, Status, Gender, Age, Ethnicity, cytopus) %>% 
                filter(Diagnosis %in% c("Control","IPF"))                
# Filter SPECTRA output to only contain cell in metadata
cell_scores <- cell_scores %>% filter(cell_id %in% rownames(metadata))
# For each patient save celltype resolution table of cell scores
lapply(unique(metadata$Sample_Name), function(current_name){
        # Get cells for current patient
        patient_cell_scores <- cell_scores %>% filter(cell_id %in% (metadata %>% filter(Sample_Name == current_name) %>% rownames()))
        # group cells by celltype
        patient_cell_scores <- patient_cell_scores %>%  left_join(metadata %>% 
                                                                        rownames_to_column(var = "cell_id") %>%
                                                                        select(cell_id, cytopus) %>% 
                                                                        unique()) %>% 
                                                        relocate(cytopus, .before=cell_id) %>% 
                                                        group_by(cytopus) %>% 
                                                        summarise_at(colnames(patient_cell_scores)[-1], mean) %>% 
                                                        rename("cell_type"="cytopus")
        write_csv(patient_cell_scores, file = paste0(output_folder, "patients_cell_scores/",current_name, ".csv"))
})
# Save patient status
patient_status <- metadata %>% select(Sample_Name, Diagnosis) %>% unique()
write_csv(patient_status, file = paste0(output_folder, "patients_status.csv"))
