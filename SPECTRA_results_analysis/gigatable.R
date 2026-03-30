# Filippo Gastaldello - 05/03/26

# Pull together the results of program activation from different methods to 
# create a consensus table

library(tidyverse)
library(jsonlite)
library(readxl)

genesets            <- c("CYTOPUS")

genesets_json_path  <- list("CYTOPUS"="/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/resources/genesets/Cytopus_unlisted_genesets.json")
program_name_column <- list("CYTOPUS"="program_short")

program_description <- read_excel("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cytopus_gene_sets/cytopus_genesets_description.xlsx", sheet = 1)

for (geneset_name in genesets) {
        
        # Load all result tables only keeping FDR and log2FC/effect size
        
        # LMM
        lmm_patient       <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/lmm/IPG_activation_lmm.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC) %>% mutate(cell_type = gsub("/","_", cell_type))
        lmm_patient_phase <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/lmm/IPG_activation_lmm_phase_corrected.csv")) %>% dplyr::filter(term == "DiagnosisIPF") %>%  dplyr::select(cell_type, program, FDR, log2FC) %>% mutate(cell_type = gsub("/","_", cell_type))
        # Welch
        welch_activation  <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/welch_test/IPG_activation_welch_test.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC) %>% mutate(cell_type = gsub("/","_", cell_type))
        welch_expression  <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/welch_test/IPG_expression_welch_test.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC) %>% mutate(cell_type = gsub("/","_", cell_type))
        # Wilcoxon-Mann-Whitney
        wmw_activation    <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/wilcoxon_test/IPG_activation_WMW_test.csv")) %>% dplyr::select(cell_type, program, FDR, effect_size) %>% mutate(cell_type = gsub("/","_", cell_type))
        wmw_expression    <- read.csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/wilcoxon_test/IPG_expression_WMW_test.csv")) %>% dplyr::select(cell_type, program, FDR, effect_size) %>% mutate(cell_type = gsub("/","_", cell_type))
        
        gigatable <- lmm_patient %>%    mutate(program_short = str_split_i(program, "-X-", 3)) %>% 
                                        relocate(program_short, .after = program) %>% 
                                        full_join(lmm_patient_phase, by = c("cell_type", "program"), suffix = c("_lmm", "_lmm_phase")) %>% 
                                        full_join(welch_activation %>% dplyr::rename("FDR_welch_activation" = "FDR", "log2FC_welch_activation" = "log2FC"), by = c("cell_type", "program")) %>% 
                                        full_join(welch_expression %>% dplyr::rename("FDR_welch_expression" = "FDR", "log2FC_welch_expression" = "log2FC", "program_short" = "program"), by = c("cell_type", "program_short")) %>% 
                                        full_join(wmw_activation %>% dplyr::rename("FDR_wmw_activation" = "FDR", "effect_size_wmw_activation" = "effect_size"), by = c("cell_type", "program")) %>% 
                                        full_join(wmw_expression %>% dplyr::rename("FDR_wmw_expression" = "FDR", "effect_size_wmw_expression" = "effect_size", "program_short" = "program"), by = c("cell_type", "program_short")) %>% 
                                        mutate(program_short = ifelse(is.na(program_short), str_split_i(program, "-X-", 3), program_short))
        
        # Load all GSEA and ORA results
        paths <- list.files(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset_name),"_genesets/manual_annotation/enrichment"),
                            pattern = "csv",
                            full.names = TRUE,
                            recursive = TRUE)
        paths <- paths[grepl(geneset_name, paths)]
        
        
        # GSEA
        
        # Read all GSEA csv, for sc, pb and both DEGs, adding padj and NES for each category
        for (DEG_set in c("sc", "pb")) {
                # Get GSEA results for DEGs according to 'DEG_set'
                if(!is_empty(paths[grepl(paste0("gsea/", DEG_set), paths)])){
                        current_paths <- paths[grepl(paste0("gsea/", DEG_set), paths)]
                        enrich_res <- lapply(current_paths, function(path){
                                # Get celltype from file name
                                celltype <- str_split_i(path, "/", -1)
                                celltype <- str_split_i(celltype, "\\.",1)
                                # Read csv, attach celltype and return
                                df <- read_csv(path)
                                return(data_frame(celltype, df$pathway, df$padj, df$ES, df$NES))
                        })
                        enrich_res <- bind_rows(enrich_res)
                }else{
                        enrich_res <- data.frame(NA, NA, NA, NA, NA)
                }
                colnames(enrich_res) <- c("cell_type", "program", paste0("padj_GSEA_", DEG_set), paste0("ES_GSEA_", DEG_set), paste0("NES_GSEA_", DEG_set))
                # Bind to gigatable
                gigatable <- gigatable %>% full_join(enrich_res %>% dplyr::rename("program_short" = "program"), by = c("cell_type", "program_short"))
        }
        
        # ORA UP
        
        # Read all ORA csv, for sc, pb and both DEGs, adding padj and NES for each category
        for (DEG_set in c("sc", "pb", "both")) {
                # Get ORA results for DEGs according to 'DEG_set'
                if(!is_empty(paths[grepl(paste0("ora/", DEG_set,".+UP"), paths)])){
                        current_paths <- paths[grepl(paste0("ora/", DEG_set,".+UP"), paths)]
                        enrich_res <- lapply(current_paths, function(path){
                                # Get celltype from file name
                                celltype <- str_split_i(path, "/", -1)
                                celltype <- str_split_i(celltype, "\\.",1)
                                celltype <- str_split_i(celltype, "_UP",1)
                                # Read csv, attach celltype and return
                                df <- read_csv(path)
                                # Actually compute gene ration
                                df <- df %>% mutate(GeneRatio = (as.numeric(str_split_i(GeneRatio, "/",1))/as.numeric(str_split_i(GeneRatio, "/",2))))
                                return(data_frame(celltype, df$ID, df$p.adjust, df$GeneRatio))
                        })
                        enrich_res <- bind_rows(enrich_res)
                }else{
                        enrich_res <- data.frame(NA, NA, NA, NA)
                }
                colnames(enrich_res) <- c("cell_type", "program", paste0("padj_ORA_UP_", DEG_set), paste0("GeneRatio_ORA_UP_", DEG_set))
                enrich_res <- enrich_res %>% dplyr::filter(.[[3]]<0.05)
                # Bind to gigatable
                gigatable <- gigatable %>% full_join(enrich_res %>% dplyr::rename("program_short" = "program"), by = c("cell_type", "program_short"))
        }
        
        # ORA DOWN
        
        # Read all ORA csv, for sc, pb and both DEGs, adding padj and NES for each category
        for (DEG_set in c("sc", "pb", "both")) {
                # Get ORA results for DEGs according to 'DEG_set'
                if(!is_empty(paths[grepl(paste0("ora/", DEG_set,".+DOWN"), paths)])){
                        current_paths <- paths[grepl(paste0("ora/", DEG_set,".+DOWN"), paths)]
                        enrich_res <- lapply(current_paths, function(path){
                                # Get celltype from file name
                                celltype <- str_split_i(path, "/", -1)
                                celltype <- str_split_i(celltype, "\\.",1)
                                celltype <- str_split_i(celltype, "_DOWN",1)
                                # Read csv, attach celltype and return
                                df <- read_csv(path)
                                # Actually compute gene ration
                                df <- df %>% mutate(GeneRatio = (as.numeric(str_split_i(GeneRatio, "/",1))/as.numeric(str_split_i(GeneRatio, "/",2))))
                                return(data_frame(celltype, df$ID, df$p.adjust, df$GeneRatio))
                        })
                        enrich_res <- bind_rows(enrich_res)
                }else{
                        enrich_res <- data.frame(NA, NA, NA, NA)
                }
                colnames(enrich_res) <- c("cell_type", "program", paste0("padj_ORA_DOWN_", DEG_set), paste0("GeneRatio_ORA_DOWN_", DEG_set))
                enrich_res <- enrich_res %>% dplyr::filter(.[[3]]<0.05)
                # Bind to gigatable
                gigatable <- gigatable %>% full_join(enrich_res %>% dplyr::rename("program_short" = "program"), by = c("cell_type", "program_short"))
        }
        
        # Add program/geneset size
        
        # Compute size of each program/set
        geneset         <- read_json(genesets_json_path[[geneset_name]])
        geneset_size    <- lapply(geneset, length)
        geneset_size_df    <- data_frame("program_short" = names(geneset_size), "geneset_size"= unname(unlist(geneset_size)))
        # Extract column containing program name for this specific iteration
        name_column  <- program_name_column[[geneset_name]]
        # Add new column
        gigatable    <- gigatable %>%
                                left_join(geneset_size_df, by = name_column) %>%
                                relocate(geneset_size, .after = any_of(name_column))
        
        # Split program name and cell type specificity
        
        
        gigatable <- gigatable %>%
                        mutate(program_specificity = str_split_i(program_short, "_",1)) %>%
                        relocate(program_specificity, .before = "program_short")
        
        # Add program topic
        gigatable <- gigatable %>% 
                        left_join(program_description %>% 
                                          select(gene_set_name, gene_set_topic),
                                  by = join_by("program_short"=="gene_set_name")) %>%
                        dplyr::rename("program_topic" = "gene_set_topic") %>% 
                        relocate(program_topic, .after = program_short)
        
        write_csv(gigatable, "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/cytopus_genesets/manual_annotation/gigatable.csv")
}