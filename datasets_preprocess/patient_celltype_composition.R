# Filippo Gastaldello - 22/01/2026

# Compute cell type composition of different patients in the dataset

library(Seurat)
library(tidyverse)
library(UpSetR)

# Load Seurat object
data <- readRDS("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/Natri_et_al.rds")
# Get metadata of patients of interest
metadata <- data@meta.data %>% 
                filter(Diagnosis == "Control" | Diagnosis == "IPF")

res <- lapply(unique(metadata$Sample_Name), function(current_sample){
        
                patient_composition <- metadata %>%
                                filter(Sample_Name==current_sample) %>%
                                select(cytopus) %>% 
                                unique()
                return(data_frame("count" = 1, "composition" = paste(sort(patient_composition$cytopus), collapse = "&")))
        })
celltype_composition <- bind_rows(res)
celltype_composition <- celltype_composition %>% 
                                group_by(composition) %>% 
                                summarize(sum(count)) %>%
                                rename("patients"="sum(count)") %>% 
                                mutate(celltype_count = str_count(composition, "&")+1)


lapply(complete_samples, function(sample){
        print(paste0(sample, ": ", length(rownames(metadata %>% filter(Sample_Name==sample))), " - ", metadata %>% filter(Sample_Name==sample) %>% select(Diagnosis) %>% unique()))
})
