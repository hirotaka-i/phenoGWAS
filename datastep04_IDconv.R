# R script to convert IDs to go along with imputed files (vcf)
# In the VCF files of following cohort, ID is "ID_ID"
# The code to check the ID format in vcf.gz file
# less /data/LNG/CORNELIS_TEMP/progression_GWAS/SCOPA/chr22.dose.vcf.gz | grep "CHROM" 

library(dplyr)
library(data.table)
for (COHORT in c("CORIELL", "DATATOP", "PARKFIT", "PICNICS", "PRECEPT", "SCOPA")){
  tmp = fread(paste("outputs/rvtest/", COHORT, ".covPC", sep = ""), sep = "\t") %>% 
    mutate(FID = paste(FID, FID, sep="_"),
           IID = paste(IID, IID, sep="_")) %>% 
    write.table(., paste("outputs/rvtest/", COHORT, ".covPC", sep = ""), sep = "\t", quote = F, row.names = F)
  print(paste(COHORT, "is converted"))
}

rm(list = ls())
print("complete")
