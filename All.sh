#!/bin/bash
# All process

# Filtering SNPs by the standard maf=0.001, Rsq03 -> output file COHORT_maf001rsq03_chr1.info/txt (rangeFile)
for dataset in "CORIELL"  "DATATOP"	"HBS" "OSLO"     "PARKFIT"  "PARKWEST" "PDBP"     "PICNICS"  "PPMI"     "PRECEPT"  "SCOPA";do
	sbatch --cpus-per-task=20 --mem=100g --time=1:00:00 SNPfilter.sh $dataset 0.001 0.3;
done
# sbach = 10


# Create PCs
##  set exclusion_region
echo '5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4' > exclusion_regions_hg19.txt
## Create PC1-5 and combine with .cov files
for dataset in "CORIELL"  "DATATOP"	"HBS" "OSLO"     "PARKFIT"  "PARKWEST" "PDBP"     "PICNICS"  "PPMI"     "PRECEPT"  "SCOPA";do
	sbatch --cpus-per-task=20 --mem=40g --time=8:00:00 PlinkPC5.sh $dataset;
done;
# sbatch = 10


# Create phenotypes files for analysis for imputed individuals
## HBS data has not genetic-phenotypic table so will be eliminated from the process
## For rvtest: .binom, .contlm .cov and MODELs.txt (outputs/rvtest/.)
## For survival: OUTCOME.COHORT.surv (outputs/survival/.)
## MODELs.txt: cohort:covs:outcomes(lm):outcomes(glm):outcomes(surv)
## PC1-5 are joined with .cov file or .surv file.
module load R
Rscript datastep01_CreatePhenoFile.R


# The IDs in imputed files are different from the original file
## CORIELL, DATATOP, PARKFIT, PARKWEST, PICNICS, PRECEPT, SCOPA -> ID_ID
## OSLO, PDBP, PPMI: OK
## HBS: HBS_PD_INVDY797MD3	HBS_PD_INVXJ077RRF


# rvtest for binomial outcomes
## Note; Oslo cohorts doesn't have the baseline characteristics. (retro/prospective combined. So doesn't have a specific baseline)
for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  model=$(echo $item | cut -d ":" -f 2)
  outcomeS=$(echo $item | cut -d ":" -f 4 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $i == "" ];then
      continue
    fi
    for chrnum in {1..22};do
      sbatch --cpus-per-task=20 --mem=40g --time=8:00:00 rvtest.sh $outcome $chrnum maf001rsq03 $dataset $model,PC1,PC2,PC3,PC4,PC5
    done
  done
done
