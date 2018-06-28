#!/bin/bash
# All process

# Create Phenofile with genetic IDs
# HBS data has not genetic-phenotypic table so will be eliminated from the process.
module load R
Rscript datastep01_CreatePhenoFile.R


# Filtering SNPs by the standard maf=0.001, Rsq03 -> output file COHORT_maf001rsq03_chr1.info/txt (rangeFile)
for dataset in "CORIELL"  "DATATOP"  "OSLO"     "PARKFIT"  "PARKWEST" "PDBP"     "PICNICS"  "PPMI"     "PRECEPT"  "SCOPA"; 
do sbatch --cpus-per-task=20 --mem=100g --time=1:00:00 SNPfilter.sh $dataset 0.001 0.3;
done



##########
# rvtest #
##########

# Devide PhenoFile for rvtest
# Outcomes file(.phoeno), ID file (.id), and covariate file (.cov) will be created.
# COVsSet.txt is a list of available COVs per cohort.
Rscript datastep03_rvtest.R

# Create PCs
##  set exclusion_region
echo '5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4' > exclusion_regions_hg19.txt
## Create PC1-5 and combine with .cov files
for item in $(tail -n +2 outputs/rvtest/COVsSet.txt);do 
	dataset=$(echo $item | cut -d ":" -f 1);
	sbatch --cpus-per-task=20 --mem=40g --time=8:00:00 PlinkPC5withCov.sh $dataset;
done;
## IDs in covPC files should be changed to the plinkr imputed data format.
### WARNING!! Resulting file names *.covPC are the same as the pre-converted file.
### So check if the file have already converted or not before running this. 
Rscript datastep04_IDconv.R



# rvtest for binomial outcomes
for item in $(tail -n +2 outputs/rvtest/COVsSet.txt);do 
	dataset=$(echo $item | cut -d ":" -f 1);
	model=$(echo $item | cut -d ":" -f 2);
	for chrnum in {1..22};do 
		sbatch --cpus-per-task=20 --mem=40g --time=8:00:00 rvtest.sh DEMENTIA $chrnum maf001rsq03 $dataset $model;
	done;
done;
