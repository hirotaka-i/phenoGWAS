#!/bin/bash
# All process

# Filtering SNPs by the standard maf=0.001, Rsq03 -> output file COHORT_maf001rsq03_chr1.info/txt (rangeFile)
for dataset in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
	sbatch --cpus-per-task=2 --mem=4g --time=6:00:00 SNPfilter.sh $dataset 0.001 0.3;
done
# sbach = 10

# More stringent filter. maf>0.05, rsq>0.8
for dataset in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
	sbatch --cpus-per-task=2 --mem=4g --time=6:00:00 SNPfilter.sh $dataset 0.05 0.8;
done


# Create PCs
##  set exclusion_region
echo '5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4' > exclusion_regions_hg19.txt
## Create PC1-5 and combine with .cov files
for dataset in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
	sbatch --cpus-per-task=2 --mem=4g --time=8:00:00 PlinkPC5.sh $dataset;
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



# rvtest for binomial outcomes
# Use swarm
# sbatch --cpus-per-task=4 --mem=8g --time=8:00:00 -> removed for swarm
rm rvtest.swarm
for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  model=$(echo $item | cut -d ":" -f 2)
  outcomeS=$(echo $item | cut -d ":" -f 4 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    for chrnum in {1..22};do
      echo  rvtest.sh $outcome $chrnum maf001rsq03 $dataset $model,PC1,PC2,PC3,PC4,PC5 >> rvtest.swarm
    done
  done
done
# 1341 jobs

# rvtest for continous outcomes
for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  model=$(echo $item | cut -d ":" -f 2)
  outcomeS=$(echo $item | cut -d ":" -f 3 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    for chrnum in {1..22};do
      echo rvtest.sh $outcome $chrnum maf001rsq03 $dataset $model,PC1,PC2,PC3,PC4,PC5 >> rvtest.swarm
    done
  done
done
# 22 jobs (PARKFIT-MMSE)

swarm -f test.swarm -g 8 -t 4 --time=1:00:00 -b 2

swarm -f test.swarm -g 8 -t 4 --time=1:00:00 -b 2 --devel

# meta-analysis
## metaanalysis per outcome
## /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/$DATASET/toMeta.GWAS.tab
rm res2meta.swarm
for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  outcomeS=$(echo $item | cut -d ":" -f 4 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    echo res2meta.sh $outcome $dataset >> res2meta.swarm
  done
done

for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  outcomeS=$(echo $item | cut -d ":" -f 3 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    echo res2meta.sh $outcome $dataset >> res2meta.swarm
  done
done

swarm -f test.swarm -g 64 -t 4  --time=8:00:00

# metal
## create metal.txt and conduct swarm
for OUTCOME in $(ls /data/LNG/Hirotaka/progGWAS/rvtest); do
 echo "
    SCHEME STDERR
    AVERAGEFREQ ON
    MINMAXFREQ ON
    " > /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt
done

for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  outcomeS=$(echo $item | cut -d ":" -f 4 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    echo "
    MARKER markerID
    ALLELE minorAllele majorAllele
    FREQ   maf
    EFFECT beta
    STDERR se
    PVALUE P
    PROCESS /data/LNG/Hirotaka/progGWAS/rvtest/$outcome/$dataset/toMeta.GWAS.tab
    " >> /data/LNG/Hirotaka/progGWAS/rvtest/$outcome/metal.txt
  done
done

for OUTCOME in $(ls /data/LNG/Hirotaka/progGWAS/rvtest); do
  echo "
  OUTFILE /data/LNG/Hirotaka/progGWAS/rvtest/$outcome/META.tbl
  ANALYZE HETEROGENEITY
  QUIT
  " >> /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt
done

rm metal.swarm
for OUTCOME in $(ls /data/LNG/Hirotaka/progGWAS/rvtest); do
	echo metal /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt >> metal.swarm
done

swarm -f metal.swarm -g 8 -t 4 --time=1:00:00 --module metal