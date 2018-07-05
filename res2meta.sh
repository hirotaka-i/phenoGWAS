#!bin/bash
# sbatch results.sh $outcome $dataset
# a combined data at /data/LNG/Hirotaka/progGWAS/rvtest/$PHENO/$DATASET/toMeta.GWAS.tab
PHENO=$1
DATASET=$2


cd /data/LNG/Hirotaka/progGWAS
rm -rf ./rvtest/$PHENO/$DATASET
mkdir -p ./rvtest/$PHENO/$DATASET
cat ./rvtest/$PHENO/chr*/$DATASET.SingleWald.assoc | grep -v 'N_INFORMATIVE' > ./rvtest/$PHENO/$DATASET/allChrs.assoc
cat ./SNPfilter/$DATASET"_"maf001rsq3_chr*.info | grep -v 'Rsq' > ./rvtest/$PHENO/$DATASET/allChrs.Info
cd ./rvtest/$PHENO/$DATASET

module load R
echo '
require(dplyr)
require(data.table)
infos <- fread(paste("allChrs.Info"))
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
assoc <- fread(paste("allChrs.assoc"))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- left_join(assoc, infos, by = c("Test" = "SNP"))
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
dat$markerID <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat$Rsq <- dat$Rsq
dat0 <- dat[,c("markerID","minorAllele","majorAllele","beta","se","maf","P", "Rsq")]
write.table(dat0, file=paste("toMeta.GWAS.tab"), quote = F, sep = "\t", row.names = F)
' > res2meta.R
Rscript --vanilla res2meta.R
