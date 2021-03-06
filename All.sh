#!/bin/bash
# All process
# set DATASETs


# Preparations
echo '
CORIELL
DATATOP
HBS
OSLO
PARKFIT
PARKWEST
PDBP
PICNICS
PPMI
PRECEPT
SCOPA
DIGPD_chip
DIGPD_neuroX' > Portfolios.txt

GNTYP_IN=/data/LNG/CORNELIS_TEMP/progression_GWAS
GNTYP_OUT=/data/LNG/Hirotaka/progGWAS

# The following filter is used
FILTER=maf001rsq3
# The following is the place for transposed SNPs cut by 20K each
FOLDER=$GNTYP_OUT/SNPfilter/$FILTER"_"20Kcut
WORKDIR=$(pwd)



mkdir outputs # outputs for phenotypes



# Create IDlist of imputed individuals
## Note delimiters of .fam files are different. Both space and tap are used. So put DATASET name at the front, take 1st column and then separate.
rm -f _IDlist_imputed.txt
for DATASET in $(cat Portfolios.txt);do
  less $GNTYP_IN/$DATASET/plink_files_hard/*.fam | sed "s/^/$DATASET:/g" | cut -f1 | cut -d " " -f1 | sed 's/:/'$'\t/g' >> _IDlist_imputed.txt
done




# QC MAP and Rsq filtering for imputerd SNPs
echo '
args <- commandArgs(trailingOnly = TRUE)
COHORT = args[1]
MAF_thres = as.numeric(args[2])
RSQ_thres = as.numeric(args[3])
GNTYP_IN = args[4]
GNTYP_OUT = arg[5]
library(data.table)
for(i in 1:22){
  print(paste(COHORT, " chr", i))
  TEXT = paste("tail -n +1", GNTYP_IN, "/", COHORT, "/chr", i, ".info.gz | gunzip", sep = "")
  data = fread(TEXT)
  dat <- subset(data, MAF >= MAF_thres & Rsq >= RSQ_thres)
  dat$chr <- tstrsplit(as.character(dat$SNP), ":")[[1]]
  dat$bp <- tstrsplit(as.character(dat$SNP), ":")[[2]]
  dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
  da <- dat[,c("SNP","ALT_Frq","Rsq")]
  da1 <- dat[,c("range")]
  LABEL = paste("_maf", substr(MAF_thres,3,nchar(MAF_thres)), "rsq", substr(RSQ_thres, 3, nchar(RSQ_thres)), "_chr", sep = "")
  write.table(da, paste(GNTYP_OUT, "/SNPfilter/", COHORT, LABEL, i, ".info",sep = ""), row.names = F, quote = F, sep = "\t")
  write.table(da1, paste(GNTYP_OUT, "/SNPfilter/", COHORT, LABEL, i, ".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}' > _SNPfilter.R

# MAF 0.001 and Rsq 0.3 or MAF 0.05 and Rsq 0.8
rm _SNPfilter.swarm
for DATASET in $(cat Portfolios.txt) ;do
  echo "
  Rscript --vanilla SNPfilter.R $DATASET 0.001 0.3 $GNTYP_IN $GNTYP_OUT
  Rscript --vanilla SNPfilter.R $DATASET 0.05 0.8 $GNTYP_IN $GNTYP_OUT" >> _SNPfilter.swarm
done

rm -rf swarm_SNPfilter
swarm -f _SNPfilter.swarm --time=0:20:00 -g 3 -p 2 -b 3 --logdir ./swarm_SNPfilter --module R




fefe
FILTER=maf001rsq3
# Transpose SNPs like the following format for the analysis in R.
## ....format....
## ID         SNP1    SNP2  ...
## Sample1    1.23    0.11
## Sample2...

# First cut each chromosome in a easy-operable number of SNPs (20K)
## "SEP" in R code will define the number. 
## Count the number of filtered SNPs in each chromosome 
rm -f _sep1.txt
for DATASET in $(cat Portfolios.txt) ;do
  for CHRNUM in {1..22};do
    echo $DATASET $CHRNUM `cat $GNTYP_OUT/SNPfilter/$DATASET"_"$FILTER"_"chr$CHRNUM.txt | wc -l` >> _sep1.txt
  done
done

## Then define the location to cut
echo '
options(scipen=999)
num = read.table("_sep1.txt")
SEP = 20000
num$V4 = num$V3/SEP
num$V5 = floor(num$V4) # The number of separation for the chromosome
t = rep(NA, nrow(num)*100) # the max snp is 120*10K so 100 would be more than enough
x = 1
for (i in 1:nrow(num)){
  for (j in 0:num$V5[i]){
    if(j==0){START=1}else{START=SEP*j} #First 2 are "ID" and "DOSE" not SNPs
    if(j==num$V5[i]){END=num$V3[i]}else{END=SEP*(j+1)-1}
    t[x] = paste(num$V1[i], ".", num$V2[i], ".", j+1, ".", START, ",", END, "p", sep = "")
    x = x + 1
  }
}
KEYS = data.frame(SEPKEY = t[!is.na(t)])
write.table(KEYS, "_sep2.txt", row.names = F, quote = F, sep = "\t")
' > _sep2.R
Rscript --vanilla _sep2.R

## Cut the SNPs with above defined locations and will use it for transposing pipeline.
mkdir $GNTYP_OUT/SNPfilter/$FILTER"_"20Kcut
FOLDER=$GNTYP_OUT/SNPfilter/$FILTER"_"20Kcut
for DATASET in $(cat Portfolios.txt) ;do
  mkdir -p $FOLDER/$DATASET
  for i in $(cat _sep2.txt | grep $DATASET);do
    CHNUM=$(echo $i | cut -d '.' -f 2)
    ITER=$(echo $i | cut -d '.' -f 3)
    AUG=$(echo $i | cut -d '.' -f 4) # location to cut. Will use for "sed -n $AUG ...."
    echo "
    cd $GNTYP_OUT/SNPfilter; \
    sed -n $AUG $GNTYP_OUT/SNPfilter/$DATASET"_"$FILTER"_"chr$CHNUM.txt | sed 's/:/\t/g' | sed 's/-/\t/g' > $FOLDER/$DATASET/chr$CHNUM.$ITER.txt.conv; \
    cd $FOLDER/$DATASET; \
    bcftools view -R chr$CHNUM.$ITER.txt.conv $GNTYP_IN/$DATASET/chr$CHNUM.dose.vcf.gz | bgzip -c > chr$CHNUM.$ITER.vcf.gz; \
    DosageConvertor --vcfDose chr$CHNUM.$ITER.vcf.gz \
      --type mach \
      --format 1 --prefix chr$CHNUM.$ITER.filter;\
    zless chr$CHNUM.$ITER.filter.mach.info | cut -f 1,2,3 | sed 's/\t/_/g' | sed '1d'> chr$CHNUM.$ITER.variant_list.txt; \
    echo DOSE > dose.txt; \
    echo ID > ID.txt; \
    cat ID.txt dose.txt chr$CHNUM.$ITER.variant_list.txt > chr$CHNUM.$ITER.final_variant_list.txt; \
    cut -f1 chr$CHNUM.$ITER.final_variant_list.txt | paste -s | sed 's/ /\t/g' > chr$CHNUM.$ITER.final_variant_list_trans.txt; \
    gunzip chr$CHNUM.$ITER.filter.mach.dose.gz; \
    cat chr$CHNUM.$ITER.filter.mach.dose | sed 's/->/\t/g' | cut -f2- > chr$CHNUM.$ITER.filter.mach.dose.format; \
    cat chr$CHNUM.$ITER.final_variant_list_trans.txt chr$CHNUM.$ITER.filter.mach.dose.format | gzip > chr$CHNUM.$ITER.trans.txt.gz; \
    rm -f chr$CHNUM.$ITER.f* chr$CHNUM.$ITER.v* chr$CHNUM.$ITER.txt* dose.txt ID.txt
    " >> _vcf_transtext.swarm
  done
done

# 5K subjobs so use 
rm -rf swarm_vcf_transtext
swarm -f _vcf_transtext.swarm -g 1 -p 2 -b 100 --time=0:04:00 --module samtools,dosageconvertor --logdir ./swarm_vcf_transtext











# Create PCs from original arrayed-data. 
##  set exclusion_region
echo '5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4' > exclusion_regions_hg19.txt

## Create PC1-5
rm _PlinkPC5.swarm
for DATASET in $(cat Portfolios.txt) ;do
  echo "
  plink --noweb --bfile $GNTYP_IN/pre_impute_data/$DATASET \
  --exclude range exclusion_regions_hg19.txt --maf 0.05 --hwe 0.00001 --geno 0.05 --indep 50 5 2 --out $GNTYP_OUT/pca/$DATASET"_"pruning; \
  plink --noweb --bfile $GNTYP_IN/pre_impute_data/$DATASET \
  --extract $GNTYP_OUT/pca/$DATASET"_"pruning.prune.in --make-bed --out $GNTYP_OUT/pca/$DATASET"_"pruned; \
  plink --noweb --bfile $GNTYP_OUT/pca/$DATASET"_"pruned --pca 5 header tabs --out $GNTYP_OUT/pca/$DATASET.pca
  " >> _PlinkPC5.swarm
done
 
rm -rf swarm_PlinkPC5
swarm -f _Plink PC5.swarm --time=0:10:00 -g 2 -p 2 -b 6 --logdir ./swarm_PlinkPC5 --module plink

# FID, IID in PC5 files are different from those in imputed files in some datasets. Correct this.
# temporary make eigenvec2 because PARKFIT returns blank if direct re-writing...
for DATASET in $(cat Portfolios.txt) ;do
  if echo "$DATASET" | grep -q '^HBS\|^OSLO\|^PDBP\|^PPMI';then
    continue
  fi
  cat $GNTYP_OUT/pca/$DATASET*.eigenvec | awk '{print $1"_"$1"\t"$2"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $GNTYP_OUT/pca/$DATASET.pca.eigenvec2
  mv $GNTYP_OUT/pca/$DATASET.pca.eigenvec2 $GNTYP_OUT/pca/$DATASET.pca.eigenvec
done



























###############################################################################################################################################
# Create phenotypes files for analysis for imputed individuals
## For rvtest (baseline): .binom, .contlm, .cov and MODELs.txt (outputs/rvtest/.)
## For rvtest (slope): .slope files and outputs/long_slope_outcomes.txt (DATASET and OUTCOME)
## For survival: OUTCOME.COHORT.surv (outputs/survival/.)
## For 
## MODELs.txt: cohort:covs:outcomes(lm):outcomes(glm):outcomes(surv)
## PC1-5 are joined with .cov file or .surv file.


# Match HBS IID and clinical ID
WORKDIR=$(pwd)
cd $GNTYP_OUT/HBS_original
plink --bfile pruned --bmerge $GNTYP_IN/pre_impute_data/HBS --make-bed --out V1
plink --bfile pruned --exclude V1-merge.missnp --make-bed --out V2
plink --bfile V2 --bmerge $GNTYP_IN/pre_impute_data/HBS --make-bed --out V3
plink --bfile V3 --genome --min 0.7
cp plink.genome $WORKDIR/data

# plink.genome is the file for ID matching
cd $WORKDIR

# Create phenotype file for rvtests and longitudinal analysis
module load R
Rscript --vanilla datastep01_CreatePhenoFile.R








# time-associated allelic effect 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4592098/ 


rm -rf outputs/long_slope
mkdir outputs/long_slope

echo '
# create data for slope analysis from longitudinal data
library(data.table)
library(dplyr)
library(lme4)

DATASETs = list.files(path = "outputs/long/") %>% tstrsplit(., "\\.") %>% .[[1]]

x = 1
AVAIL_OUTCOMEs = data.frame(DatasetOutcome=NA) # record OUTCOMEs for later analysis

for (DATASET in DATASETs){
  # Read data
  cohort=fread(paste("outputs/long/", DATASET, ".cont", sep = "")) %>% arrange(IID, YEARfDIAG) %>% mutate(ID = paste(FID, IID, sep="_"))

  # OUTCOME list, COVs list
  temp = fread("outputs/MODELs.txt", sep=":", skip=1)
  OUTCOMEs = temp %>% filter(V1==DATASET) %>% with(tstrsplit(V6, ",")) %>% unlist

  # Keep the time varying variables for regression
  ## the treatment can slow the disease progression
  COVs_trans = names(cohort)[grep("^YEARfDIAG|^DOPA|^AGONIST", names(cohort))]


  # OUTCOME iteratioin
  ## Prepare OUTPUT file
  OUTPUT = cohort %>% select(FID, IID, FATID, MATID, SEX) %>% mutate(ID = paste(FID, IID, sep="_")) %>% distinct(IID, .keep_all = T)
  for (OUTCOME in OUTCOMEs){
    cohort$OUTCOME = cohort[, OUTCOME]
    cohort_temp = cohort %>% filter(!is.na(OUTCOME))

    # Transform the data to the orthogonal to the cross-sectional space
    cohort_trans = cohort_temp %>% select(c("ID", "OUTCOME", COVs_trans)) %>% arrange(ID, YEARfDIAG)
    IDs = unique(cohort_trans$ID)
    trans.func = function(i){
      xi = cohort_trans %>% filter(ID == IDs[i]) %>% select(-ID)
      xi = as.matrix(xi)
      if(nrow(xi) >1){
        A = cumsum(rep(1, nrow(xi)))
        A1 = poly(A, degree = length(A)-1)
        transxi = t(A1) %*% xi
        transxi[abs(transxi) < 0.0000000001] = 0 # put 0
        transxi = data.frame(ID = IDs[i], transxi)
        return(transxi)
      }else{
        return(rep(NA, length(c("OUTCOME", COVs_trans))+1))
      }
    }

    temp = lapply(1:length(IDs), trans.func)
    transdata = do.call(rbind, temp)

    # Determine slopes for individuals, note we do not need intercept (-1 |ID)
    if (length(COVs_trans)>1){             
      testMod = try(lmer(paste("OUTCOME ~ ", paste(setdiff(COVs_trans, "YEARfDIAG"), collapse="+"), "- 1 + (YEARfDIAG -1 | ID)"), data = transdata), silent=T)
    }else{
      testMod = try(lmer(OUTCOME ~ 1 + (YEARfDIAG -1 | ID), data = transdata), silent=T)
    }
    if(class(testMod)[1]!="try-error"){
      blups = ranef(testMod)$ID
      names(blups)=paste(OUTCOME, "slope", sep="_")
      blups$ID = rownames(blups)
      OUTPUT = left_join(OUTPUT, blups, by = "ID")
      AVAIL_OUTCOMEs[x, 1] = paste(DATASET, paste(OUTCOME, "slope", sep="_"), sep = ":")
      x = x + 1
      print(paste(DATASET, OUTCOME, "DONE"))
    }
  }
  write.table(OUTPUT, paste("outputs/long_slope", "/", DATASET, ".slope", sep=""), row.names = F, quote = F, sep = "\t")
}
write.table(AVAIL_OUTCOMEs, "outputs/long_slope_outcomes.txt", row.names = F, quote = F, sep = "\t")
' > _slope.R


Rscript --vanilla _slope.R
# the R.code create slope file for rvtest as well as available outcome files.


























##################################################################################################################################################
# Analysis (cohort level)
# cross-sectional -> rvtests
# longitudinal -> R


###############################################
# rvtest - for outcomes per individual ########
###############################################
rm -f _rvtest.swarm _rvtest_done.txt
# For binomial outcomes (logistic regression)
for item in $(tail -n +2 outputs/MODELs.txt);do 
  DATASET=$(echo $item | cut -d ":" -f 1)
  MODEL=$(echo $item | cut -d ":" -f 2)
  OUTCOMEs=$(echo $item | cut -d ":" -f 4 |tr , "\t")
  for OUTCOME in $OUTCOMEs;do
    if [ $OUTCOME == "" ];then
      continue
    fi
    echo "rvtest_binom,$DATASET,$OUTCOME" >> _rvtest_done.txt
    for CHNUM in {1..22};do
      echo "
      mkdir -p $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM;\
      rvtest --noweb --hide-covar --rangeFile $GNTYP_OUT/SNPfilter/$DATASET"_"$FILTER"_"chr$CHNUM.txt \
      --inVcf $GNTYP_IN/$DATASET/chr$CHNUM.dose.vcf.gz \
      --pheno outputs/rvtest/$DATASET.binom --pheno-name $OUTCOME \
      --covar outputs/rvtest/$DATASET.cov --covar-name $MODEL \
      --out $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM/$DATASET --single wald
      ">>_rvtest.swarm
    done
  done
done

# For  continous outcome (linear regression)
for item in $(tail -n +2 outputs/MODELs.txt);do 
  DATASET=$(echo $item | cut -d ":" -f 1)
  MODEL=$(echo $item | cut -d ":" -f 2)
  OUTCOMEs=$(echo $item | cut -d ":" -f 3 |tr , "\t") #Only here changed
  for OUTCOME in $OUTCOMEs;do
    if [ $OUTCOME == "" ];then
      continue
    fi
    echo "rvtest_cont,$DATASET,$OUTCOME" >> _rvtest_done.txt
    for CHNUM in {1..22};do
      echo "
      mkdir -p $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM;\
      rvtest --noweb --hide-covar --rangeFile $GNTYP_OUT/SNPfilter/$DATASET"_"$FILTER"_"chr$CHNUM.txt \
      --inVcf $GNTYP_IN/$DATASET/chr$CHNUM.dose.vcf.gz \
      --pheno outputs/rvtest/$DATASET.binom --pheno-name $OUTCOME \
      --covar outputs/rvtest/$DATASET.cov --covar-name $MODEL \
      --out $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM/$DATASET --single wald
      ">>_rvtest.swarm
    done
  done
done

# For slope, (linear regression without covariates)
for item in $(tail -n +2 outputs/long_slope_outcomes.txt);do
  OUTCOME=$(echo $item | cut -d ":" -f 2)
  DATASET=$(echo $item | cut -d ":" -f 1)
  echo "rvtest_slope,$DATASET,$OUTCOME" >> _rvtest_done.txt
  for CHNUM in {1..22};do
    echo "
    mkdir -p $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM;\
    rvtest --noweb --hide-covar --rangeFile $GNTYP_OUT/SNPfilter/$DATASET"_"$FILTER"_"chr$CHNUM.txt \
  --inVcf $GNTYP_IN/$DATASET/chr$CHNUM.dose.vcf.gz \
  --pheno outputs/long_slope/$DATASET.slope --pheno-name $OUTCOME \
  --out $GNTYP_OUT/rvtest/$OUTCOME/chr$CHNUM/$DATASET --single wald
    " >> _rvtest.swarm
  done
done


rm -rf swarm_rvtest
swarm -f _rvtest.swarm --time=6:00:00 -p 2 -b 8 --logdir ./swarm_rvtest --module rvtests --devel
# maf001rsq3.. longetst time was 4:30:00 (b4) so 2hr for each job would be enough

# Oslo didn't finish in 2:00
# rm -rf swarm_rvtest
# swarm -f ../OSLOimcomplete.txt --time=6:00:00 -p 2 --logdir ./swarm_rvtest --module rvtests --devel






###########################################
# Longitudinal outcomes ###################
# Binomial -> survival analysis ###########
# Continous -> linear mixed model  ########
###########################################
# Cox analysis
# Code for analysis
# # e.g. Rscript --vanilla cox_single.R CONST DATATOP 22 9  $GNTYP_OUT/SNPfilter/$FILTER_20Kcut
echo '
# cox analysis
args <- commandArgs(trailingOnly = TRUE)
# set variables
OUTCOME = args[1]
DATASET = args[2]
chrNUM =  args[3]
ITER =  args[4]
FOLDER = args[5] 

# Cox model for GWAS
library(data.table)
library(dplyr)
library(survival)

# Read data
## DATASET with TSTART and TSTOP
cohort=fread(paste("outputs/surv/", OUTCOME, ".", DATASET, ".surv", sep = "")) %>% arrange(IID, TSTART)
COVs = paste(names(cohort)[-c(1:3,ncol(cohort)-5:7)], collapse=" + ") # dropb c(ID, TSTART, OUTCOME, BEGIN, END, TSTOP)

## Imputed data
SNPset = fread(paste("zcat -f ", FOLDER, "/", DATASET, "/", chrNUM, ".", ITER, ".trans.txt.gz", sep=""))
SNPs = names(SNPset)[-(1:2)] # 1 ID, 2 DOSE, SNP name starts from 3
## Merge
cohort_snp = left_join(cohort, SNPset, by = c("IID"="ID"))
cohort_snp$SurvObj1 = with(cohort_snp, Surv(TSTART, TSTOP, OUTCOME == 1))

surv.listfunc = function(x){
  # Models
  MODEL = paste("SurvObj1 ~", "`", SNPs[x], "`+", COVs, sep = "")
  testCox = try(coxph(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testCox)[1]=="try-error"){
    sumstat=c(SNPs[x], "NoConverge", rep(NA,4))
  }else{
    temp= summary(testCox)$coefficients
    if(grep(SNPs[x], rownames(temp)) %>% length == 0){ # In this case, SNP is dropeed from the model
      sumstat=c(SNPs[x], "NoVforSNP", rep(NA, 4))
    }else{
      RES = temp[1,]
      EVENT_OBS = paste(testCox$nevent, testCox$n, sep="_")
      sumstat <- c(SNPs[x], EVENT_OBS, as.numeric(RES[4]), RES[1], RES[3], RES[5])
    }
  }
  return(sumstat)
}

temp = lapply(1:length(SNPs), surv.listfunc)
temp2 = do.call(rbind, temp)
attributes(temp2)$dimnames[[2]]=c("SNP", "EVENT_OBS", "TEST", "Beta", "SE", "Pvalue")
NEWDIR = paste("/data/LNG/Hirotaka/progGWAS/surv/", OUTCOME, "/", chrNUM, sep = "")
dir.create(NEWDIR, recursive = T, showWarnings = F)
write.table(temp2, paste(NEWDIR, "/", DATASET, ".", ITER, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
print("finish")
' > _surv.R


FOLDER=$GNTYP_OUT/SNPfilter/$FILTER"_"20Kcut
rm _surv.swarm _surv_done.txt
for i in $(ls outputs/surv/);do 
  OUTCOME=$(echo $i | cut -d '.' -f 1)
  DATASET=$(echo $i | cut -d '.' -f 2)
  echo "surv,$DATASET,$OUTCOME" >> _surv_done.txt
  for j in $(ls $FOLDER/$DATASET);do
    chrNUM=$(echo $j | cut -d '.' -f 1)
    ITER=$(echo $j | cut -d '.' -f 2)
    echo "Rscript --vanilla _surv.R $OUTCOME $DATASET $chrNUM $ITER $FOLDER">> _surv.swarm
  done
done

rm -rf swarm_surv
swarm -f _surv.swarm -g 1.5 -p 2 -b 100 --time=0:18:00 --module R --logdir ./swarm_surv
# The initial setting of -b 130 got timeouts aroung 102nd jobs at 24:00:00


# Retrieve the imcomplete jobs and re-submit
JOBNUM=6270580
rm _fail_$JOBNUM*
jobhist $JOBNUM | grep 'TIMEOUT' | cut -d " " -f1 | cut -d "_" -f2 > _fail_$JOBNUM.txt
for i in $(cat _fail_$JOBNUM.txt);do
  cat /spin1/swarm/$USER/$JOBNUM/cmd.$i"_"* >> _fail"_"$JOBNUM"_"temp.txt # when -p option on, two jobs per cpu. e.g. cmd.11_0/1 
done
# Take the input and only retireve the one without outputs. (Only for the program not updating until the job completes.)
while read item;do
  OUTCOME=$(echo $item | cut -d " " -f 4)
  DATASET=$(echo $item | cut -d " " -f 5)
  chrNUM=$(echo $item | cut -d " " -f 6)
  ITER=$(echo $item | cut -d " " -f 7)
  if ! test -e "$GNTYP_OUT/surv/$OUTCOME/$chrNUM/$DATASET.$ITER.txt";then
    echo "$item" >> _fail"_"$JOBNUM.swarm
  fi
done < _fail"_"$JOBNUM"_"temp.txt
swarm -f _fail"_"$JOBNUM.swarm -g 1.5 -p 2 -b 3 --time=0:18:00 --module R --logdir ./swarm_surv --devel


















# GLMM
# Code for analysis
# # e.g. Rscript --vanilla cox_single.R HY.DATATOP.22.9.txt $GNTYP_OUT/SNPfilter/CONV/maf001rsq3_10Kcut
echo '
# longitudinal analysis (for cross-sectional section)
args <- commandArgs(trailingOnly = TRUE)
OUTCOME = args[1]
DATASET = args[2]
chrNUM =  args[3]
ITER =  args[4]
FOLDER = args[5] 

# Cox model for GWAS
library(data.table)
library(dplyr)
library(lme4)

# Read data
## DATASET with TSTART and TSTOP
cohort=fread(paste("outputs/long/", DATASET, ".cont", sep = "")) %>% arrange(IID, YEARfDIAG) %>% mutate(ID = paste(FID, IID, sep="_"))
COVs = names(cohort)[grep("^FEMALE|^YEARSEDUC|^FAMILY_HISTORY|^AAO|^YEARfDIAG|^DOPA|^AGNOIST|^PC", names(cohort))]

## Imputed data
SNPset = fread(paste("zcat -f ", FOLDER, "/", DATASET, "/", chrNUM, ".", ITER, ".trans.txt.gz", sep=""))
SNPs = names(SNPset)[-(1:2)] # 1 ID, 2 DOSE, SNP name starts from 3
## Merge
cohort_snp = left_join(cohort, SNPset, by = c("IID"="ID")) # !!!!!!!NO PCs!!!!!!!ID IID IS IT RIGHT!!!!!!!!!!!!!!????????????


glmm.listfunc = function(x){
  # Models
  MODEL = paste(OUTCOME, "~", "`", SNPs[x], "`+", paste(COVs, collapse="+"), "+(1|ID)", sep = "")
  testLmer = try(lmer(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testLmer)[1]=="try-error"){
    sumstat=rep("DROP",6)
  }else{
    temp = summary(testLmer)
    temp1 = temp$coefficients
    if(grep(SNPs[x], rownames(temp1)) %>% length == 0){ # In this case, SNP is dropeed from the model
      sumstat=rep(NA,6)
    }else{
      RES = temp1[2,] # The first row is intercept
      PV_APPROX = 2 * pnorm(abs(RES[3]), lower.tail=F) # df is enough large for approximation
      OBS_N = paste(length(temp$residuals), "_", temp$ngrps, sep="")
      sumstat <- c(SNPs[x], OBS_N, RES[3], RES[1], RES[2], PV_APPROX)
    }
  }
  return(sumstat)
}

temp = lapply(1:length(SNPs), glmm.listfunc)
temp2 = do.call(rbind, temp)

attributes(temp2)$dimnames[[2]]=c("SNP", "OBS_N", "Tvalue", "Beta", "SE", "Pv_approx")
NEWDIR = paste("/data/LNG/Hirotaka/progGWAS/long/", OUTCOME, "/", chrNUM, sep = "")
dir.create(NEWDIR, recursive = T, showWarnings = F)
write.table(temp2, paste(NEWDIR, "/", DATASET, ".", ITER, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
' > _long.R

# Analysis
rm _long.swarm _long_done.txt
for I in $(tail -n +2 outputs/MODELs.txt);do
  OUTCOMEs=$(echo $I | cut -d ":" -f 6 | tr , "\t")
  DATASET=$(echo $I | cut -d ":" -f 1)
  for OUTCOME in $OUTCOMEs; do
    echo "long,$DATASET,$OUTCOME" >> _long_done.txt
    for J in $(ls $FOLDER/$DATASET/);do
      chrNUM=$(echo $J | cut -d "." -f1)
      ITER=$(echo $J | cut -d "." -f2)
      echo "Rscript --vanilla _long.R $OUTCOME $DATASET $chrNUM $ITER $FOLDER" >> _long.swarm
    done
  done
done

swarm -f _long.swarm -g 1.5 -p 2 -b 24 --time=2:00:00 --module R --logdir ./swarm_long --devel

# For glmm_slope
head -n 28 _long.R > _long_slope.R
echo '
  MODEL = paste(OUTCOME, "~", "`", SNPs[x], "`+ `", SNPs[x], "`*YEARfDIAG + ", paste(COVs, collapse="+"), "+(1 + YEARfDIAG|ID)", sep = "")
  lmer(eval(parse(text = MODEL)), data = cohort_snp)
  
  
  testLmer = try(lmer(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testLmer)[1]=="try-error"){
    sumstat=rep("DROP",6)
  }else{
    temp = summary(testLmer)
    temp1 = temp$coefficients
    GET = nrow(temp1)
    if(grep(SNPs[x], rownames(temp1)[GET]) %>% length == 0){ # In this case, SNP is dropeed from the model
      sumstat=rep(NA,6)
    }else{
      RES = temp1[GET,] # The first row is intercept
      PV_APPROX = 2 * pnorm(abs(RES[3]), lower.tail=F) # df is enough large for approximation
      OBS_N = paste(length(temp$residuals), "_", temp$ngrps, sep="")
      sumstat <- c(SNPs[x], OBS_N, RES[3], RES[1], RES[2], PV_APPROX)
    }
  }
  return(sumstat)
}

temp = lapply(1:length(SNPs), glmm.listfunc)
temp2 = do.call(rbind, temp)

attributes(temp2)$dimnames[[2]]=c("SNP", "OBS_N", "Tvalue", "Beta", "SE", "Pv_approx")
NEWDIR = paste("/data/LNG/Hirotaka/progGWAS/long/", OUTCOME, "_slope/", chrNUM, sep = "")
dir.create(NEWDIR, recursive = T, showWarnings = F)
write.table(temp2, paste(NEWDIR, "/", DATASET, ".", ITER, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
' >> _long_slope.R

rm _long_slope.swarm _long_slope_done.txt
for I in $(tail -n +2 outputs/MODELs.txt);do
  OUTCOMEs=$(echo $I | cut -d ":" -f 6 | tr , "\t")
  DATASET=$(echo $I | cut -d ":" -f 1)
  for OUTCOME in $OUTCOMEs; do
    echo "long,$DATASET,$OUTCOME" >> _long_slope_done.txt
    for J in $(ls $FOLDER/$DATASET/);do
      chrNUM=$(echo $J | cut -d "." -f1)
      ITER=$(echo $J | cut -d "." -f2)
      echo "Rscript --vanilla _long_slope.R $OUTCOME $DATASET $chrNUM $ITER $FOLDER" >> _long_slope.swarm
    done
  done
done



swarm -f _long_slope.swarm -g 3 -p 2 -b 24 --time=3:00:00 --module R --logdir ./swarm_long_slope --devel



























##################################################################################################################################################
# meta-analysis
## meta_analysis per outcome
# first prepare the output files



# Reference Info file
## Input file, chr*.info.gz is not always referencing the major allele
## To reference major allele, create allChrs.Info and use ALT_FREQ in the later step
for DATASET in $(cat Portfolios.txt) ;do
  cat $GNTYP_OUT/SNPfilter/"$DATASET"_"$FILTER"_chr*.info | grep -v 'Rsq' > $GNTYP_OUT/SNPfilter/"$DATASET"_"$FILTER"_allChrs.Info
done



## For METAL, format the analysis results like the followings.
# markerID  minorAllele majorAllele beta  se  maf P Rsq
# chr10:64972 A G 0.201144  1.52208 0.00942 0.894865  0.55998
# chr10:65030 A C 0.0774858 0.27058 0.44703 0.774595  0.40803
echo '
args <- commandArgs(trailingOnly = TRUE)
OUTCOME=args[1]
DATASET=args[2]
FILTER=paste("_", args[3], "_", sep = "") # ="_maf001rsq3_"

require(dplyr)
require(data.table)
infos <- fread(paste("../SNPfilter", "/", DATASET, FILTER, "allChrs.Info", sep=""))
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
assoc <- fread(paste(OUTCOME, "/", DATASET, "/allChrs.assoc", sep=""))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- left_join(assoc, infos, by = c("Test" = "SNP")) %>% 
  filter(Beta < 5 & Beta > -5) %>% 
  filter(!is.na(Pvalue)) %>% 
  mutate(chr = paste("chr",CHROM, sep = ""),
         ALT_Frq=as.numeric(ALT_Frq),
         Beta=as.numeric(Beta)) %>% 
  mutate(markerID = paste(chr,POS, sep = ":"),
         noflip = if_else(ALT_Frq <= 0.5, TRUE, FALSE),
         minorAllele = if_else(noflip, as.character(ALT), as.character(REF)),
         majorAllele = if_else(noflip, as.character(REF), as.character(ALT)),
         beta = ifelse(noflip, Beta, Beta*(-1)),
         se = SE,
         P = Pvalue,
         maf = if_else(noflip, ALT_Frq, 1 - ALT_Frq)) %>% 
  select(markerID, minorAllele, majorAllele, beta, se, maf, P, Rsq)
write.table(data, file=paste(OUTCOME, "/", DATASET, "/toMeta.GWAS.tab", sep = ""), quote = F, sep = "\t", row.names = F)
' > _metal_input.R

# Grip association test results from rvtest and convert it for the metal input
WORKDIR=$(pwd)
rm _metal_input.swarm
for I in $(tail -n +1 _rvtest_done.txt) ;do
  DATASET=$(echo $I | cut -d "," -f 2)
  OUTCOME=$(echo $I | cut -d "," -f 3)
  echo "
  cd $GNTYP_OUT/rvtest;\
  rm -rf $OUTCOME/$DATASET;\
  mkdir  $OUTCOME/$DATASET;\
  cat $OUTCOME/chr*/$DATASET.SingleWald.assoc | grep -v N_INFORMATIVE > $OUTCOME/$DATASET/allChrs.assoc;\
  Rscript --vanilla $WORKDIR/_metal_input.R $OUTCOME $DATASET $FILTER
  " >> _metal_input.swarm
done


# Do the same thing for cox results. 
for I in $(tail -n +1 _surv_done.txt) ;do
  DATASET=$(echo $I | cut -d "," -f 2)
  OUTCOME=$(echo $I | cut -d "," -f 3)
  echo "
  cd $GNTYP_OUT/surv;\
  rm -rf $OUTCOME/$DATASET;\
  mkdir  $OUTCOME/$DATASET;\
  cat $OUTCOME/chr*/$DATASET*.txt | grep -v Pvalue | sed 's/:/\t/g' | sed 's/_/\t/' | sed 's/_/\t/' |\
  awk '"'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$1":"$2"\t"$7"\t"$8"\t"$9}'"' > $OUTCOME/$DATASET/allChrs.assoc;\
  Rscript --vanilla $WORKDIR/_metal_input.R $OUTCOME $DATASET $FILTER
  " >> _metal_input.swarm
done


# Also for glmm resutls
for I in $(tail -n +1 _long_done.txt) ;do
  DATASET=$(echo $I | cut -d "," -f 2)
  OUTCOME=$(echo $I | cut -d "," -f 3)
  echo "
  cd $GNTYP_OUT/long;\
  rm -rf $OUTCOME/$DATASET;\
  mkdir  $OUTCOME/$DATASET;\
  cat $OUTCOME/chr*/$DATASET*.txt | grep -v Pvalue | sed 's/:/\t/g' | sed 's/_/\t/' | sed 's/_/\t/' |\
  awk '"'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$1":"$2"\t"$7"\t"$8"\t"$9}'"' > $OUTCOME/$DATASET/allChrs.assoc;\
  Rscript --vanilla $WORKDIR/_metal_input.R $OUTCOME $DATASET $FILTER
  " >> _metal_input.swarm
done
swarm -f _metal_input.swarm -g 10 -p 2 -b 60 --time=0:05:00 --module R --logdir ./swarm_metal --devel


# Create metal.txt 
rm -rf $GNTYP_OUT/meta 
mkdir $GNTYP_OUT/meta

###################LONG_slope, LONG, count number!! ####################
for ANALYSIS in "rvtest" "long" "surv"; do
  for OUTCOME in $(ls $GNTYP_OUT/$ANALYSIS); do
    mkdir -p $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME
    echo "
    SCHEME STDERR
    AVERAGEFREQ ON
    MINMAXFREQ ON
    " > $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/metal.txt
    for DATASET in $(tail -n +1 _"$ANALYSIS"_done.txt | grep "$OUTCOME"$ | cut -d ',' -f2);do
      echo "
      MARKER markerID
      ALLELE minorAllele majorAllele
      FREQ   maf
      EFFECT beta
      STDERR se
      PVALUE P
      PROCESS $GNTYP_OUT/$ANALYSIS/$OUTCOME/$DATASET/toMeta.GWAS.tab
      " >> $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/metal.txt
      echo "$DATASET" >> $GNTYP_OUT/meta/temp.txt
    done
    echo "
    OUTFILE $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/meta .tbl
    ANALYZE HETEROGENEITY
    QUIT
    " >> $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/metal.txt
    DATASETS=$(paste -s $GNTYP_OUT/meta/temp.txt | sed 's/\t/,/g')
    echo "$ANALYSIS"___"$OUTCOME"___"$DATASETS" >>  _meta_done.txt
    rm $GNTYP_OUT/meta/temp.txt
  done
done

# Note MMSE for PARKFIT was only cross-sectional, so need to combine this with glmm model analyses
# Create the updated version and rename originals.
cp $GNTYP_OUT/meta/long/MMSE/metal.txt $GNTYP_OUT/meta/long/MMSE/metal_original.txt
head -n 41  $GNTYP_OUT/meta/long/MMSE/metal_original.txt > $GNTYP_OUT/meta/long/MMSE/metal.txt
sed -n '7,13p' $GNTYP_OUT/meta/rvtest/MMSE/metal.txt >> $GNTYP_OUT/meta/long/MMSE/metal.txt
tail -n +42  $GNTYP_OUT/meta/long/MMSE/metal_original.txt >> $GNTYP_OUT/meta/long/MMSE/metal.txt
mv $GNTYP_OUT/meta/rvtest/MMSE/metal.txt $GNTYP_OUT/meta/rvtest/MMSE/metal_1cohort.txt 


rm _metal.swarm
for FILE in $(ls $GNTYP_OUT/meta/*/*/metal.txt); do
  echo "metal $FILE" >> _metal.swarm
done

# conduct swarm
swarm -f _metal.swarm -g 6 -p 2 -b 2 --time=1:00:00 --module metal --logdir ./swarm_metal --module metal --devel


head -n +2 _metal.swarm > test.swarm
swarm -f test.swarm -g 6 -p 2 --time=1:00:00 --module metal --logdir ./swarm_metal --module metal --devel





# Filter the meta-analysis results
# The Df threshold for the meta-analysis. using SNPs available more than 50% of the cohorts.
## Note the calculation doesn't have decimals.
rm -f _metathres.txt
for i in $(cat _meta_done.txt);do
  N_COHORTs=$(echo $i | sed 's/___/:/g' | cut -d ":" -f 3 | sed 's/,/\n/g' | wc -l)
  echo $[(N_COHORTs) / 2 ]  >> _metathres.txt
done
paste -d '-' _meta_done.txt _metathres.txt > _meta_done_thres.txt
## Filtering HetIsq <80 HetDf > half of the number of the available cohorts.
ANALYSIS=surv
for i in $(grep ^"$ANALYSIS" _meta_done_thres.txt);do
  OUTCOME=$(echo $i | sed 's/___/:/g' | cut -d ":" -f 2)
  THRES=$(echo $i | cut -d "-" -f 2)
  awk '$14>'$THRES' {print $0}' $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/meta1.tbl | awk '$12<80' | sort -gk 10 > $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/meta2.tbl
done


ANALYSIS=rvtest
for i in $(grep ^"$ANALYSIS" _meta_done_thres.txt);do
  OUTCOME=$(echo $i | sed 's/___/:/g' | cut -d ":" -f 2)
  THRES=$(echo $i | cut -d "-" -f 2)
  awk '$14>'$THRES' {print $0}' $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/meta1.tbl | awk '$12<80' | sort -gk 10 > $GNTYP_OUT/meta/$ANALYSIS/$OUTCOME/meta2.tbl
done



# DATA CHECKING COMMAND
sort -gk 7 /Volumes/Hirotaka/progGWAS/rvtest/MMSE_slope/PPMI/toMeta.GWAS.tab > sorted_test.tab
grep 'chr1:' sorted_test.tab | head -n 20


grep '155205634\|155206037\|155206167' sorted_test2.tab