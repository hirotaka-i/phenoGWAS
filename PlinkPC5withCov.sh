#!/bin/bash
DATASET=$1

module load plink
plink --noweb --keep outputs/rvtest/$DATASET.id --bfile /data/LNG/CORNELIS_TEMP/progression_GWAS/pre_impute_data/$DATASET \
--exclude range exclusion_regions_hg19.txt --maf 0.05 --hwe 0.00001 --geno 0.05 --indep 50 5 2 --out /data/LNG/Hirotaka/progGWAS/$DATASET"_"pruning
plink --noweb --keep outputs/rvtest/$DATASET.id --bfile /data/LNG/CORNELIS_TEMP/progression_GWAS/pre_impute_data/$DATASET \
--extract /data/LNG/Hirotaka/progGWAS/$DATASET"_"pruning.prune.in --make-bed --out /data/LNG/Hirotaka/progGWAS/$DATASET"_"pruned
plink --noweb --keep outputs/rvtest/$DATASET.id --bfile /data/LNG/Hirotaka/progGWAS/$DATASET"_"pruned --pca 5 header tabs --out /data/LNG/Hirotaka/progGWAS/$DATASET"_"pca
### Join other covariates file with PC1-5
tail -n +1 /data/LNG/Hirotaka/progGWAS/$DATASET"_"pca.eigenvec | awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$7}' > /data/LNG/Hirotaka/progGWAS/$DATASET"_"pca.noID
paste -d '\t' outputs/rvtest/$DATASET.cov /data/LNG/Hirotaka/progGWAS/$DATASET"_"pca.noID > outputs/rvtest/$DATASET.covPC

