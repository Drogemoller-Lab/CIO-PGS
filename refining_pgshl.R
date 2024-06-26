###Author: Deanne Nixie Miao
###Date: April 21, 2024
###Goal: Filter variants mapping to genes within the six filters

#####Installilng packages
install.packages("babelgene")
install.packages("janitor")
library(babelgene)
library(tidyverse)
library(dplyr)
library(janitor)

###Reading datasets
#Read in relevant data
kalra_sumstats_GRCh38 <- fread("/home/projects/archive/cisplatin/PGS/kalra_sumstats_vep_annotation_GRCh38.tsv", head =TRUE)
kalra_hg38_vep <- fread("/home/projects/archive/cisplatin/VEP/kalra_hg38_for_VEP.aa_results", header = TRUE)
mouse_milo_adjusted <- fread("/home/projects/cio_pgs/new_sc_gene_lists/DESeq2_per_cluster_merged_pval005_treat_up_dn_milo_3logFC.txt")

###Clean datasets
kalra_hg38_vep <- clean_names(kalra_hg38_vep)
kalra_hg38_vep <- kalra_hg38_vep[,c(1,4)]
kalra_hg38_vep$number_uploaded_variation <- gsub("chr", "", kalra_hg38_vep$number_uploaded_variation)
kalra_hg38_vep$number_uploaded_variation <- gsub("_", ":", kalra_hg38_vep$number_uploaded_variation)
kalra_hg38_vep$number_uploaded_variation <- gsub("/", ":", kalra_hg38_vep$number_uploaded_variation)

kalra_sumstats_GRCh38$id <- paste(kalra_sumstats_GRCh38$chromosome, kalra_sumstats_GRCh38$position,
                                  kalra_sumstats_GRCh38$effect_allele, kalra_sumstats_GRCh38$non_effect_allele,
                                  sep = ":")
kalra_sumstats_GRCh38$id <- gsub("chr", "", kalra_sumstats_GRCh38$id)

mouse_milo_adjusted <- mouse_milo_adjusted[,-1]
sum(duplicated(mouse_milo_adjusted$gene)) #23
mouse_milo_adjusted <- mouse_milo_adjusted %>% distinct(gene, .keep_all = TRUE) #159 left

### Filter
mouse_milo_adjusted_human_data <- orthologs(genes = mouse_milo_adjusted$gene, 
                                            species = "Mus musculus", human = F)

mouse_milo_adjusted_human_data <- mouse_milo_adjusted_human_data[,c("human_ensembl", "human_symbol","symbol")]
mouse_milo_adjusted_human_data <- unique(mouse_milo_adjusted_human_data)
mouse_milo_adjusted_human_data <- na.omit(mouse_milo_adjusted_human_data)

filtered_kalra_deseq_milo_adjusted <- kalra_hg38_vep[kalra_hg38_vep$gene %in% 
                                                       mouse_milo_adjusted_human_data$human_ensembl, ]
filtered_kalra_deseq_milo_adjusted <- filtered_kalra_deseq_milo_adjusted %>% distinct(number_uploaded_variation, .keep_all = TRUE)
filtered_kalra_deseq_milo_adjusted <- merge(filtered_kalra_deseq_milo_adjusted, kalra_sumstats_GRCh38,
                                            by.x = "number_uploaded_variation",
                                            by.y = "id",
                                            all.x = TRUE)

###Make sure that there are no ambiguous variants 
filtered_kalra_deseq_milo_adjusted <- filtered_kalra_deseq_milo_adjusted %>%
  filter(!((effect_allele == "T" & non_effect_allele == "A") |
             (effect_allele == "A" & non_effect_allele == "T") |
             (effect_allele == "G" & non_effect_allele == "C") |
             (effect_allele == "C" & non_effect_allele == "G")))

##### Save data
write.table(filtered_kalra_deseq_milo_adjusted, "/home/projects/archive/cisplatin/PGS/results_v4/filtered_kalra_deseq_milo_adjusted.txt", quote = FALSE, sep = "\t", row.names = FALSE)


