#!/bin/bash
# rvtest for slope. Not that it doesn't have covariates now. 
#e.g. rvtest3.sh DEMENTIA 22 maf001rsq03 PPMI
PHENO=$1
CHNUM=$2
THRES=$3
DATASET=$4 
mkdir -p /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/chr$CHNUM
rvtest --noweb --hide-covar --rangeFile /data/LNG/Hirotaka/progGWAS/SNPfilter/$DATASET"_"$THRES"_"chr$CHNUM.txt \
--inVcf /data/LNG/CORNELIS_TEMP/progression_GWAS/$DATASET/chr$CHNUM.dose.vcf.gz \
--pheno outputs/long_slope/$DATASET.slope --pheno-name $PHENO \
--out /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/chr$CHNUM/$DATASET --single wald