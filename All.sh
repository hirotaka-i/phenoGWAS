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
swarm -f SNPfilter.swarm --time=0:20:00 -g 3 -p 2 -b 3 --logdir ./swarm_SNPfilter --module R --devel


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
  OUTFILE /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME .tbl
  ANALYZE HETEROGENEITY
  QUIT
  " >> /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt
done

rm metal.swarm
for OUTCOME in $(ls /data/LNG/Hirotaka/progGWAS/rvtest); do
	echo metal /data/LNG/Hirotaka/progGWAS/rvtest/$OUTCOME/metal.txt >> metal.swarm
done

swarm -f test.swarm -g 6 -p 2 -b 2 --time=1:00:00 --module metal --logdir ./swarm_metal --module metal

####################################################################################################

# Survival analysis 
## Need to use R for using Cox model with time-varying covariates 
### Convert SNPs file to 
## format:
# ID    SNP1    SNP2
# Sample1    1.23    0.11
# etc


DATASET=$1
CHNUM=$2

#sbatch --cpus-per-task=20 --mem=100g --time=12:00:00 convert_HRC_impute_to_mach_dosage_NEW.sh PPMI 22

module load samtools
module load dosageconvertor
# sinteractive --mem=100g --cpus-per-task=20
cd /data/LNG/CORNELIS_TEMP/progression_GWAS/$DATASET

# filter variants for R2 and MAF
cp   all_filtered_sites.txt
sed -i ‘s/:/\t/g’ all_filtered_sites.txt
sed -i ‘s/-/\t/g’ all_filtered_sites.txt

# output file is .mach.dose.gz
# info file output is remove.mach.info
bcftools view -R all_filtered_sites.txt chr$CHNUM.dose.vcf.gz | bgzip -c > $CHNUM.vcf.gz
DosageConvertor --vcfDose $CHNUM.vcf.gz \
    --type mach \
    --format 1 --prefix $CHNUM.filter

# Next make sure to handle with duplicate chr:bp variant
# meaning merge with A1 and A2 and prep to add header
zless $CHNUM.filter.mach.info | cut -f 1,2,3 > variant_list.txt
sed -i ‘s/\t/_/g’ variant_list.txt
sed -i ‘1d’ variant_list.txt
echo DOSE > dose.txt
echo ID > ID.txt
cat ID.txt dose.txt variant_list.txt > final_variant_list.txt

# transpose
awk '
{
   for (i=1; i<=NF; i++)  {
       a[NR,i] = $i
   }
}
NF>p { p = NF }
END {    
   for(j=1; j<=p; j++) {
       str=a[1,j]
       for(i=2; i<=NR; i++){
           str=str” “a[i,j];
       }
       print str
   }
}' final_variant_list.txt > final_variant_list_trans.txt
sed -i ‘s/ /\t/g’ final_variant_list_trans.txt

# next add header to dosage file
gunzip $CHNUM.filter.mach.dose.gz
sed -i ‘s/->/\t/g’ $CHNUM.filter.mach.dose
cut -f2- $CHNUM.filter.mach.dose > output 
cat final_variant_list_trans.txt output | gzip > $CHNUM.filter.dosage_format.txt.gz
rm all_filtered_sites.txt
rm $CHNUM.filter.mach.info
rm $CHNUM.vcf.gz
rm $CHNUM.filter.mach.dose
rm variant_list.txt
rm dose.txt
rm ID.txt
rm final_variant_list.txt
rm final_variant_list_trans.txt
rm output