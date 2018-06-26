#!/bin/bash
# rvtest for binomial outcomes
# e.g. sbatch --cpus-per-task=20 --mem=100g --time=2:00:00 rvtest.sh PPMI 22 DEMENTIA
PHENO=$1
CHNUM=$2
THRES=$3
DATASET=$4 
MODEL=$5

module load rvtests
rvtest --noweb --hide-covar --rangeFile /data/LNG/Hirotaka/progGWAS/$DATASET"_"$THRES"_"chr$CHNUM.txt \
--inVcf /data/LNG/CORNELIS_TEMP/progression_GWAS/$DATASET/chr$CHNUM.dose.vcf.gz \
--pheno outputs/rvtest/$DATASET.pheno --pheno-name $PHENO \
--covar outputs/rvtest/$DATASET.covPC --covar-name $MODEL \
--out /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO"_"$DATASET.chr$CHNUM --single wald,score