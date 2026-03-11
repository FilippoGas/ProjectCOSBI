# Filippo Gastaldello - 05/03/26

# Pull together the results of program activation from different methods to 
# create a consensus table

library(tidyverse)

geneset <- "CYTOPUS"

# Load all result tables only keeping FDR and log2FC/effect size

# LMM
lmm_patient       <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/lmm/IPG_activation_lmm.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC)
lmm_patient_phase <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/lmm/IPG_activation_lmm_phase_corrected.csv")) %>% dplyr::filter(term == "DiagnosisIPF") %>%  dplyr::select(cell_type, program, FDR, log2FC)
# Welch
welch_activation  <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/welch_test/IPG_activation_welch_test.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC)
welch_expression  <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/welch_test/IPG_expression_welch_test.csv")) %>% dplyr::select(cell_type, program, FDR, log2FC)
# Wilcoxon-Mann-Whitney
wmw_activation    <- read_csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/wilcoxon_test/IPG_activation_WMW_test.csv")) %>% dplyr::select(cell_type, program, FDR, effect_size)
wmw_expression    <- read.csv(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/wilcoxon_test/IPG_expression_WMW_test.csv")) %>% dplyr::select(cell_type, program, FDR, effect_size)

gigatable <- lmm_patient %>%    mutate(program_short = str_split_i(program, "-X-", 3)) %>% 
                                relocate(program_short, .after = program) %>% 
                                full_join(lmm_patient_phase, by = c("cell_type", "program"), suffix = c("_lmm", "_lmm_phase")) %>% 
                                full_join(welch_activation %>% dplyr::rename("FDR_welch_activation" = "FDR", "log2FC_welch_activation" = "log2FC"), by = c("cell_type", "program")) %>% 
                                full_join(welch_expression %>% dplyr::rename("FDR_welch_expression" = "FDR", "log2FC_welch_expression" = "log2FC", "program_short" = "program"), by = c("cell_type", "program_short")) %>% 
                                full_join(wmw_activation %>% dplyr::rename("FDR_wmw_activation" = "FDR", "effect_size_wmw_activation" = "effect_size"), by = c("cell_type", "program")) %>% 
                                full_join(wmw_expression %>% dplyr::rename("FDR_wmw_expression" = "FDR", "effect_size_wmw_expression" = "effect_size", "program_short" = "program"), by = c("cell_type", "program_short"))

# Load all GSEA and ORA results
paths <- list.files(paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/",tolower(geneset),"_genesets/manual_annotation/enrichment"),
                    pattern = "csv",
                    full.names = TRUE,
                    recursive = TRUE)
paths <- paths[grepl(geneset, paths)]


# GSEA

# Read all GSEA csv, for sc, pb and both DEGs, adding padj and NES for each category
for (DEG_set in c("sc", "pb", "both")) {
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
        colnames(enrich_res) <- c("cell_type", "program", paste0("padj_GWAS_", DEG_set), paste0("ES_GWAS_", DEG_set), paste0("NES_GWAS_", DEG_set))
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
        colnames(enrich_res) <- c("cell_type", "program", paste0("padj_ORA_UP_", DEG_set), paste0("GeneRatio_ORA_up_", DEG_set))
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
        colnames(enrich_res) <- c("cell_type", "program", paste0("padj_ORA_DOWN_", DEG_set), paste0("GeneRatio_ORA_DOWN_", DEG_set))
        enrich_res <- enrich_res %>% dplyr::filter(.[[3]]<0.05)
        # Bind to gigatable
        gigatable <- gigatable %>% full_join(enrich_res %>% dplyr::rename("program_short" = "program"), by = c("cell_type", "program_short"))
}

write_csv(gigatable, "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/cytopus_genesets/manual_annotation/gigatable.csv")
