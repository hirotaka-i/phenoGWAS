#!/bin/bash
DATASET=$1

module load plink
plink --noweb --bfile /data/LNG/CORNELIS_TEMP/progression_GWAS/pre_impute_data/$DATASET \
--exclude range exclusion_regions_hg19.txt --maf 0.05 --hwe 0.00001 --geno 0.05 --indep 50 5 2 --out /data/LNG/Hirotaka/progGWAS/pca/$DATASET"_"pruning
plink --noweb --bfile /data/LNG/CORNELIS_TEMP/progression_GWAS/pre_impute_data/$DATASET \
--extract /data/LNG/Hirotaka/progGWAS/pca/$DATASET"_"pruning.prune.in --make-bed --out /data/LNG/Hirotaka/progGWAS/pca/$DATASET"_"pruned
plink --noweb --bfile /data/LNG/Hirotaka/progGWAS/pca/$DATASET"_"pruned --pca 5 header tabs --out /data/LNG/Hirotaka/progGWAS/pca/$DATASET.pca