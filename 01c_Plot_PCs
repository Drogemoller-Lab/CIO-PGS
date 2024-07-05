#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

#load data
cio_2024_pca <- read.table("~/Desktop/Xu_LDpruned_pca.eigenvec", header = F)
colnames(cio_2024_pca) <- c("sample","sample2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#plot and save data
ggplot(cio_2024_pca, aes(x=PC1, y=PC2)) + geom_point()
ggsave("CIO_pc1_pc2.tiff", 
       height=5.36, width=7.58, units='in', dpi=300)

#create a column specifying if the data is European or not
cio_2024_pca_ancestry <- cio_2024_pca %>%
  mutate(ancestry=case_when(PC1>=0 & PC2>=-0.005 ~ "Eu",
                            TRUE ~ "Other"))

#plot data with European ancestry indivdiauls in a different colour
ggplot(cio_2024_pca_ancestry, aes(x=PC1, y=PC2)) + geom_point(aes(colour=ancestry))
ggsave("CIO_pc1_pc2_by_ancestry.tiff", 
       height=5.36, width=7.58, units='in', dpi=300)
#subset Eu only data
cio_2024_Eu <- subset(cio_2024_pca_ancestry, ancestry=="Eu")

