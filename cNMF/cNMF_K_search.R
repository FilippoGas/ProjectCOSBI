# Filippo Gastaldello - 02/07/26
#
# Run cNMF on a scRNAseq dataset

library(Matrix)
library(Seurat)
library(tidyverse)

###############
#    PATHS    #  
###############

sc_dataset_path  <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/Natri_et_al_clean.rds"
mtx_folder       <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cNMF/mtx/"
prepare_folder   <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cNMF/prepare/"
factorize_folder <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/cNMF/factorize/"

##############
#    MAIN    #
##############


# load seurat object
sc_dataset <- read_rds(sc_dataset_path)

# Extract count matrix
counts     <- sc_dataset@assays$RNA$counts
barcodes   <- colnames(counts)
gene_names <- rownames(counts)


# Save count matrix, barcodes and gene names as a 10x-Genomics-formatted mtx directory
if (!dir.exists(mtx_folder)){dir.create(mtx_folder, recursive = TRUE)}
writeMM(counts, paste0(mtx_folder, "matrix.mtx"))
write.table(as.data.frame(barcodes), paste0(mtx_folder, "barcodes.tsv"), col.names = FALSE, row.names = FALSE, sep = "\t")
features <- data.frame("gene_id"=gene_names, "gene_names"=gene_names, type="Gene Expression")
write.table(as.data.frame(features), paste0(mtx_folder, "genes.tsv"), col.names = FALSE, row.names = FALSE, sep = "\t")

# Run the prepare step of cNMF

run_name <- "Natri_et_al"
if (!dir.exists(prepare_folder)){dir.create(prepare_folder, recursive = TRUE)}
cmd = paste("cnmf prepare",
            "   --output-dir", prepare_folder,
            "   --name", run_name,
            "   -c", paste0(mtx_folder, "matrix.mtx"),
            "   --max-nmf-iter 1000",
            "   -k", paste(seq(20,210,10), collapse=" "),
            "   --n-iter 20",
            sep = " ")
system(cmd)

# Run the factorization step
n_workers = 30
if (!dir.exists(factorize_folder)){dir.create(factorize_folder, recursive = TRUE)}
cmd = paste("nohup cnmf factorize",
            "   --output-dir", factorize_folder,
            "   --name", run_name,
            "   --worker-index", seq(0,n_workers-1,1),
            "   --total-workers", n_workers,
            "&",
            sep=" ")
