#!/bin/bash
# All process

# Filtering SNPs by the standard maf=0.001, Rsq0.3 -> output file COHORT_maf001rsq3_chr1.info/txt (rangeFile)
# Filtering SNPs by the standard maf=0.05, Rsq0.8 -> output file COHORT_maf001rsq8_chr1.info/txt (rangeFile)
rm SNPfilter.swarm
for dataset in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
  echo "Rscript --vanilla SNPfilter.R $dataset 0.001 0.3
  Rscript --vanilla SNPfilter.R $dataset 0.05 0.8" >> SNPfilter.swarm
done

rm -rf swarm_SNPfilter
swarm -f SNPfilter.swarm --time=0:20:00 -g 3 -p 2 -b 3 --logdir ./swarm_SNPfilter --module R


# Create PCs
##  set exclusion_region
echo '5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4' > exclusion_regions_hg19.txt

## Create PC1-5 and combine with .cov files
rm PlinkPC5.swarm
for dataset in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
  echo "bash PlinkPC5.sh $dataset" >> PlinkPC5.swarm
done
 
rm -rf swarm_PlinkPC5
swarm -f PlinkPC5.swarm --time=1:00:00 -g 2 -p 2 -b 6 --logdir ./swarm_PlinkPC5 --module plink

# Create phenotypes files for analysis for imputed individuals
## HBS data has not genetic-phenotypic table so will be eliminated from the process
## For rvtest: .binom, .contlm .cov and MODELs.txt (outputs/rvtest/.)
## For survival: OUTCOME.COHORT.surv (outputs/survival/.)
## MODELs.txt: cohort:covs:outcomes(lm):outcomes(glm):outcomes(surv)
## PC1-5 are joined with .cov file or .surv file.
module load R
Rscript --vanilla datastep01_CreatePhenoFile.R

# rvtest for binomial outcomes
# Use swarm
# sbatch --cpus-per-task=4 --mem=16g --time=8:00:00 -> removed for swarm
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
      echo  "bash rvtest.sh $outcome $chrnum maf001rsq3 $dataset $model,PC1,PC2,PC3,PC4,PC5" >> rvtest.swarm
    done
  done
done
# 1341 jobs

# rvtest for continous outcomes
## maf0.001 rsq 0.3 
for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  model=$(echo $item | cut -d ":" -f 2)
  outcomeS=$(echo $item | cut -d ":" -f 3 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    for chrnum in {1..22};do
      echo "bash rvtest2.sh $outcome $chrnum maf001rsq3 $dataset $model,PC1,PC2,PC3,PC4,PC5" >> rvtest.swarm
    done
  done
done
# 22 jobs (PARKFIT-MMSE)

rm -rf swarm_rvtest
swarm -f rvtest.swarm --time=2:00:00 -p 2 -b 4 --logdir ./swarm_rvtest
# maf001rsq3.. longetst time was 4:30:00 (b4) so 2hr for each job would be enough



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
    echo "bash res2meta.sh $outcome $dataset" >> res2meta.swarm
  done
done

for item in $(tail -n +2 outputs/MODELs.txt);do 
  dataset=$(echo $item | cut -d ":" -f 1)
  outcomeS=$(echo $item | cut -d ":" -f 3 |tr , "\t")
  for outcome in $outcomeS;do
    if [ $outcome == "" ];then
      continue
    fi
    echo "bash res2meta.sh $outcome $dataset" >> res2meta.swarm
  done
done

rm -rf swarm_res2meta
swarm -f res2meta.swarm -g 8 -p 2 -b 10 --time=2:00:00 --logdir ./swarm_res2meta
# 30 min. can reduce number in -b. 

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
  OUTFILE /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/meta .tbl
  ANALYZE HETEROGENEITY
  QUIT
  " >> /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt
done

rm metal.swarm
for OUTCOME in $(ls /data/LNG/Hirotaka/progGWAS/rvtest); do
  echo metal /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt >> metal.swarm
done

# MMSE Only one cohort. Cannot do metal at this moment.
grep -v "MMSE" metal.swarm > metal.swarm2

# conduct swarm
swarm -f test.swarm -g 6 -p 2 -b 2 --time=1:00:00 --module metal --logdir ./swarm_metal --module metal
####################################################################################################
# Survival analysis 
## Need to use R for using Cox model with time-varying covariates 
### Convert SNPs file to 
## format:
# ID    SNP1    SNP2
# Sample1    1.23    0.11
# etc


# Only retain the last file
rm vcf_transtext.swarm
rm -rf /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/$FILTER
FILTER=maf001rsq3
for DATASET in "CORIELL" "DATATOP" "HBS" "OSLO" "PARKFIT" "PARKWEST" "PDBP" "PICNICS" "PPMI" "PRECEPT" "SCOPA";do
  mkdir -p /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/$FILTER/$DATASET;
  for CHNUM in {1..22};do
    echo " \
    cd /data/LNG/Hirotaka/progGWAS/SNPfilter; \
    cat $DATASET"_"$FILTER"_"chr$CHNUM.txt | sed 's/:/\t/g' | sed 's/-/\t/g' > CONV/$FILTER/$DATASET/chr$CHNUM.txt.conv; \
    cd CONV/$FILTER/$DATASET; \
    bcftools view -R chr$CHNUM.txt.conv /data/LNG/CORNELIS_TEMP/progression_GWAS/$DATASET/chr$CHNUM.dose.vcf.gz | bgzip -c > chr$CHNUM.vcf.gz; \
    DosageConvertor --vcfDose chr$CHNUM.vcf.gz \
      --type mach \
      --format 1 --prefix chr$CHNUM.filter;\
    zless chr$CHNUM.filter.mach.info | cut -f 1,2,3 | sed 's/\t/_/g' | sed '1d'> chr$CHNUM.variant_list.txt; \
    echo DOSE > dose.txt; \
    echo ID > ID.txt; \
    cat ID.txt dose.txt chr$CHNUM.variant_list.txt > chr$CHNUM.final_variant_list.txt; \
    cut -f1 chr$CHNUM.final_variant_list.txt | paste -s | sed 's/ /\t/g' > chr$CHNUM.final_variant_list_trans.txt; \
    gunzip chr$CHNUM.filter.mach.dose.gz; \
    cat chr$CHNUM.filter.mach.dose | sed 's/->/\t/g' | cut -f2- > chr$CHNUM.filter.mach.dose.format; \
    cat chr$CHNUM.final_variant_list_trans.txt chr$CHNUM.filter.mach.dose.format | gzip > chr$CHNUM.trans.txt.gz; \
    rm -f chr$CHNUM.f* chr$CHNUM.v* chr$CHNUM.txt* dose.txt ID.txt
    " >> vcf_transtext.swarm
  done
done

rm -rf swarm_vcf_transtext
swarm -f vcf_transtext.swarm -g 3 -p 2 -b 4  --time=3:00:00 --module samtools,dosageconvertor --logdir ./swarm_vcf_transtext

# Code for analysis
echo '
# cox analysis
args <- commandArgs(trailingOnly = TRUE)
FILENAME = args[1]
CHRNUM = args[2]
FILTER = args[3]
N_CORE = as.integer(args[4])-1

# # e.g. Rscript --vanilla cox.R DEMENTIA.SCOPA.surv 22 maf001rsq3 2
# FILENAME = "DEMENTIA.CORIELL.surv"
# CHRNUM = 22
# FILTER = "maf001rsq3"
# N_CORE = 2

# set variables
OUTCOME = strsplit(FILENAME, "\\.")[[1]][1]
DATASET = strsplit(FILENAME, "\\.")[[1]][2]

# Cox model for GWAS
library(data.table)
library(dplyr)
library(survival)
library(parallel)

# Read data
## DATASET with TSTART and TSTOP
cohort=fread(paste("outputs/surv", FILENAME, sep = "/")) %>% arrange(ID, TSTART)
COVs = paste(names(cohort)[-c(1:3,ncol(cohort)-5:7)], collapse=" + ") # dropb c(ID, TSTART, OUTCOME, BEGIN, END, TSTOP)
print(COVs)
## Imputed data
READ_LOCATION = paste("zcat -f /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/", FILTER, "/", DATASET, "/chr", CHRNUM, ".trans.txt.gz", sep = "")
SNPset = fread(READ_LOCATION)
SNPs = names(SNPset)[-(1:2)] # 1 ID, 2 DOSE, SNP name starts from 3
## Merge
cohort_snp = left_join(cohort, SNPset, by = "ID")
cohort_snp$SurvObj1 = with(cohort_snp, Surv(TSTART, TSTOP, OUTCOME == 1))

surv.listfunc = function(x){
  # Models
  MODEL = paste("SurvObj1 ~", "`", SNPs[x], "`+", COVs, sep = "")
  testCox = try(coxph(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testCox)[1]=="try-error"){
    sumstat=rep(NA,6)
  }else{
    RES = summary(testCox)$coefficients[1,]
    EVENT_OBS = paste(testCox$nevent, testCox$n, sep="_")
    sumstat <- c(SNPs[x], EVENT_OBS, as.numeric(RES[4]), RES[1], RES[3], RES[5])
  }
  return(sumstat)
}

temp = mclapply(1:length(SNPs), surv.listfunc, mc.cores=N_CORE)
temp2 = do.call(rbind, temp)

# temp = mclapply(1:length(SNPs), surv.listfunc, mc.cores=N_CORE)
# temp2 = do.call(rbind, temp)
attributes(temp2)$dimnames[[2]]=c("SNP", "EVENT_OBS", "TEST", "Beta", "SE", "Pvalue")
NEWDIR = paste("/data/LNG/Hirotaka/progGWAS/surv/", OUTCOME, "/chr", CHRNUM, sep = "")
dir.create(NEWDIR, recursive = T, showWarnings = F)
write.table(temp2, paste(NEWDIR, "/", DATASET, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
' > cox_mclapply.R

# Analysis
rm -rf /data/LNG/Hirotaka/progGWAS/surv
rm -f cox.swarm

# FILTER=maf001rsq3
# N_CORE only works with 1 at this point ?
for i in $(ls outputs/surv/);do
  for j in {1..22};do
    echo Rscript --vanilla cox.R $i $j $FILTER 8>> cox.swarm
  done
done

swarm -f cox.swarm -g 50 -t 8 --time=4:00:00 --module R --logdir ./swarm_cox






# Different stragtegy. Just analyze with the single core lapply per 100K SNPs ##################################
# Separate chr$CHRNUM.trans.txt.gz
## get the number of separation
#### Number of SNPs at each allele
rm -f sep1.txt
for COHORT in $(ls  /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/maf001rsq3);do
  for CHRNUM in {1..22};do
    echo $COHORT $CHRNUM `less /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/maf001rsq3/$COHORT/chr$CHRNUM.trans.txt.gz | head -1 | tr '\t' '\n' | wc -l` >> sep1.txt
  done
done

echo '
options(scipen=999)
num = read.table("sep1.txt")
SEP = 10000
num$V4 = num$V3/SEP
num$V5 = floor(num$V4) # Will iterate with approximately 10K SNPs
t = rep(NA, nrow(num)*200) # the max snp is 120*10K
x = 1
for (i in 1:nrow(num)){
  for (j in 0:num$V5[i]){
    if(j==0){START=1}else{START=paste("1,","2,",SEP*j, sep="")} #First 2 are "ID" and "DOSE" not SNPs
    if(num$V5[i]==j){END=num$V3[i]}else{END=SEP*(j+1)-1}
    t[x] = paste(num$V1[i], ".", num$V2[i], ".", j+1, ".", START, "-", END, sep = "")
    x = x + 1
  }
}
KEYS = data.frame(SEPKEY = t[!is.na(t)])
write.table(KEYS, "sep2.txt", row.names = F, quote = F, sep = "\t")
' > sep2.R
Rscript --vanilla sep2.R

# FILTER=maf001rsq3
mkdir /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/$FILTER"_"10Kcut
for i in $(tail -n +2 sep2.txt);do
  DATASET=$(echo $i | cut -d '.' -f 1)
  CHRNUM=$(echo $i | cut -d '.' -f 2)
  ITER=$(echo $i | cut -d '.' -f 3)
  AUG=$(echo $i | cut -d '.' -f 4)
  less /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/$FILTER/$DATASET/chr$CHRNUM.trans.txt.gz | cut -f $AUG >> /data/LNG/Hirotaka/progGWAS/SNPfilter/CONV/$FILTER"_"10Kcut/$DATASET.$CHRNUM.$ITER.txt
done

