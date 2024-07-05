#!/bin/bash
PCA=/home/drogemob@med.umanitoba.ca/deanne/pca

GCTA=/opt/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static 
PLINK=/opt/plink-1.09/plink
DATA=/home/projects/archive/cisplatin/PCA/PCA_v2

mkdir -p $PCA

$PLINK --bfile $DATA/Xu_all_chr --maf 0.01 \
	--allow-extra-chr --allow-no-sex -indep-pairwise 100 20 0.25 \
	--out $PCA/Xu_LD

$PLINK --bfile $DATA/Xu_all_chr \
	--allow-extra-chr --allow-no-sex --extract $PCA/Xu_LD.prune.in \
	--make-bed --out $PCA/Xu_LDpruned

# Perform PCA 
$GCTA --bfile $PCA/Xu_LDpruned \
        --thread-num 16 --make-grm \
        --out $PCA/Xu_LDpruned_kinship

$GCTA --grm $PCA/Xu_LDpruned_kinship \
        --thread-num 16 --pca 10 \
        --out $PCA/Xu_LDpruned_pca

