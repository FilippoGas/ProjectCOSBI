# Filippo Gastaldello - 11/06/2026
#
# Count the number of unknown programs in the spectra results at different values of lambda 

library(tidyverse)

input_dir        <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA/full/manual_annotation/"
res <- lapply(
        list.files(input_dir, pattern = "cell_scores", full.names = TRUE),
        function(df_path){
                df <- read_csv(df_path)
                count <- as.numeric(length(grep("\\d+$", colnames(df)[2:length(colnames(df))], value = TRUE)))
                return(data.frame("lam"=df_path, "count"=count))
        })
unknown_programs_count <- bind_rows(res)

# Manually modify lam column
unknown_programs_count$lam <- c(0.001, 0.01, 0.5, 0.1)


p <- unknown_programs_count %>% 
        ggplot(aes(x=as.factor(lam), y =count)) +
                geom_bar(
                        stat  = "identity",
                        fill  = "orange",
                        alpha = 0.7,
                        width = 0.3
                        ) +
                theme_minimal() +
                theme(
                        axis.text.y      = element_blank(),
                        panel.grid.minor = element_blank()
                        ) +
                geom_text(
                        aes(label=count),
                        vjust = -1,
                        color = "black",
                        size  = 4
                        )+
                labs(
                        title = "Number of unknown programs in the cell score matrix at different values of the lambda parameter"
                        ) +
                xlab(
                        "Lambda"
                        ) +
                ylab(
                        "Unknown programs"
                )
ggsave(plot     = p,
       filename = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/COSBI/data/IPF_datasets/Natri_et_al/SPECTRA_results_analysis/full/cytopus_genesets/manual_annotation/unknown_programs/unknown_programs_count_vs_lambda.pdf",
       device   = "pdf",
       width    = 9,
       height   = 6)
