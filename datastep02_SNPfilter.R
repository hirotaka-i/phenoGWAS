# First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  cat("Supply COHORTS NAME")
}else{
#  eval(parse(text=args[1]))
  COHORT = args[1]
  MAF_thres = as.numeric(args[2])
  RSQ_thres = as.numeric(args[3])
  library(data.table)
  for(i in 1:22){
    print(paste(COHORT, " chr", i))
    TEXT = paste("fread('tail -n +1 /data/LNG/CORNELIS_TEMP/progression_GWAS/", COHORT, "/chr", i, ".info.gz | gunzip')", sep = "")
    data = eval(parse(text = TEXT))
    dat <- subset(data, MAF >= MAF_thres & Rsq >= RSQ_thres)
    dat$chr <- tstrsplit(as.character(dat$SNP), ":")[[1]]
    dat$bp <- tstrsplit(as.character(dat$SNP), ":")[[2]]
    dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
    da <- dat[,c("SNP","ALT_Frq","Rsq")]
    da1 <- dat[,c("range")]
    LABEL = paste("_maf", substr(MAF_thres,3,nchar(MAF_thres)), "rsq0", substr(RSQ_thres, 3, nchar(RSQ_thres)), "_chr", sep = "")
    write.table(da, paste("/data/LNG/Hirotaka/progGWAS/SNPfilter", COHORT, LABEL, i, ".info",sep = ""), row.names = F, quote = F, sep = "\t")
    write.table(da1, paste("/data/LNG/Hirotaka/progGWAS/SNPfilter", COHORT, LABEL, i, ".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
  }
}


rm(list=ls())