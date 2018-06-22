# First read in the arguments listed at the command line
args=(commandArgs(TRUE))

if(length(args)==0){
  cat("Supply COHORTS NAME LIKE-> '--args COHORT=\\\"COHORT_NAME\\\"'")
}else{
  eval(parse(text=args[[1]]))
  library(data.table)
  for(i in 1:22){
    print(paste(COHORT, " chr", i))
    TEXT = paste("fread('tail -n +1 /data/LNG/CORNELIS_TEMP/progression_GWAS/", COHORT, "/chr", i, ".info.gz | gunzip')", sep = "")
    data = eval(parse(text = TEXT))
    dat <- subset(data, MAF >= 0.001 & Rsq >= 0.30)
    dat$chr <- tstrsplit(as.character(dat$SNP), ":")[[1]]
    dat$bp <- tstrsplit(as.character(dat$SNP), ":")[[2]]
    dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
    da <- dat[,c("SNP","ALT_Frq","Rsq")]
    da1 <- dat[,c("range")]
    write.table(da, paste("outputs/rvtest/", COHORT, "_maf001rsq03_chr",i,".info",sep = ""), row.names = F, quote = F, sep = "\t")
    write.table(da1, paste("outputs/rvtest/", COHORT, "_maf001rsq03_chr",i,".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
  }
}


