# datastep2
# Create files for rvtest (logistic regression)

# .phenofiles
library(data.table)
library(dplyr)
# Read data
ltg_all=fread('outputs/PhenoFile.csv') %>% data.frame
names(ltg_all)
BinomT = names(ltg_all)[7:17]
BinomT
COVs = names(ltg_all)[c(5:6, 32, 34:37)]
COVs

COHORTs = unique(ltg_all$COHORT)
COVs = c('FEMALE', 'YEARSEDUC', "FAMILY_HISTORY", 'AAO', "BLDfDIAG",  'DOPA', 'AGONIST')
COVs_set = data.frame(COHORT_MODEL=rep(NA, length(COHORTs)))
dir.create("outputs/rvtest", showWarnings = F)
for(i in 1:length(COHORTs)){
  cohort = ltg_all %>% filter(TSTART ==0 & COHORT == COHORTs[i]) %>% 
    arrange(ID) %>% 
    mutate(FID = ID,
           IID = ID,
           FATID = 0,
           MATID = 0,
           SEX = 1+FEMALE) %>% # 0 missing, 1 male, 2 female
    mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1))) # 0 missing, 1 control, 2 case
  IDs = cohort %>% select(FID, IID)
  pheno = cohort %>% select(FID, IID, FATID, MATID, SEX, BinomT)
  cov = cohort %>% select(FID, IID, FATID, MATID, SEX, COVs) # Cov can use 0 not as missing.
  # filter COVARIATES with all the same values or all NAs
  temp = cov %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  COVs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA or Var=0
  cov_new = cov %>% select(FID, IID, FATID, MATID, COVs_i)
  COVs_i_woSex = setdiff(COVs_i, "SEX")
  write.table(pheno, paste("outputs/rvtest/", COHORTs[i], ".pheno", sep = ""), row.names = F, quote = F)
  write.table(cov_new, paste("outputs/rvtest/", COHORTs[i], ".cov", sep = ""), row.names = F, quote = F)
  write.table(IDs, paste("outputs/rvtest/", COHORTs[i], ".id", sep = ""), row.names = F, quote = F)
  COVs_set[i,] = paste(COHORTs[i], paste(COVs_i_woSex, collapse = ","), sep=":")
}
write.table(COVs_set, "outputs/rvtest/COVsSet.txt", row.names = F, quote = F)
print("complete")

rm(list =ls())
