#!/bin/bash


# e.g. 
# for dataset in "CORIELL"  "DATATOP"  "OSLO"     "PARKFIT"  "PARKWEST" "PDBP"     "PICNICS"  "PPMI"     "PRECEPT"  "SCOPA"; 
# do sbatch --cpus-per-task=20 --mem=100g --time=1:00:00 SNPfilter.sh $dataset;
# done

DATASET=$1 
MAF=$2
RSQ=$3
module load R
Rscript datastep02_SNPfilter.R $DATASET $MAF $RSQ
