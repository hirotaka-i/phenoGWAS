#!/bin/bash
# rvtest for binomial outcomes or single measured continous outcomes 
# e.g. sbatch --cpus-per-task=20 --mem=100g --time=2:00:00 rvtest.sh DEMENTIA 22 maf001rsq03 PPMI FEMALE,YEARSEDUC,FAMILY_HISTORY,AAO,BLDfDIAG,(PC1,PC2,PC3..)
PHENO=$1
CHNUM=$2
THRES=$3
DATASET=$4 
MODEL=$5


mkdir -p /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/chr$CHNUM

module load rvtests
rvtest --noweb --hide-covar --rangeFile /data/LNG/Hirotaka/progGWAS/SNPfilter/$DATASET"_"$THRES"_"chr$CHNUM.txt \
--inVcf /data/LNG/CORNELIS_TEMP/progression_GWAS/$DATASET/chr$CHNUM.dose.vcf.gz \
--pheno outputs/rvtest/$DATASET.binom --pheno-name $PHENO \
--covar outputs/rvtest/$DATASET.cov --covar-name $MODEL \
--out /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/chr$CHNUM/$DATASET --single wald,score
