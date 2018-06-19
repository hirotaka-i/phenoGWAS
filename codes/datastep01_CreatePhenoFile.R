# Phenofile
# Keep IDs available in both ltg_all and imputed files.

setwd("..")
# install.packages("openxlsx")
require(openxlsx)
require(dplyr)
require(data.table)
require(zoo)
###################################
# Add phenotype ID column for the imputed individuals
# Sheet 1 all except oslo
# Sheet 2 Oslo
###################################
temp1 = read.xlsx("data/progression_GWAS_imputed_IDs.xlsx", sheet = 1) %>% 
  arrange(DATASET, ID)
temp2 = temp1 %>% 
  group_by(DATASET, ID) %>% 
  mutate(PID = case_when(
    DATASET != "HBS"~ strsplit(ID, "_")[[1]][1],
    strsplit(ID, "_")[[1]][2] == "PD" ~ strsplit(ID, "_")[[1]][3],
    TRUE~ strsplit(ID, "_")[[1]][2]
  )) %>% data.frame
temp3 = read.xlsx("data/progression_GWAS_imputed_IDs.xlsx", sheet = 2) %>% 
  arrange(DATASET, ID) %>% 
  rename(PID = Clincal_ID)
temp4 = bind_rows(temp2, temp3) %>% 
  arrange(DATASET, ID) %>% 
  rename(COHORT = DATASET)
write.csv(temp4, "ID_imputed.csv", row.names = F)

#############################
# PhenoFile with genotyping data
# PhenoFile should be in PDcohorts/codes. Before retrieve, make sure the file is latest
#############################
temp1 = fread("ID_imputed.csv")
Pheno = fread(paste("../PDcohorts/codes/","Pheno13_", Sys.Date(), ".csv", sep = ""), stringsAsFactors = F) %>% 
  rename(COHORT = STUDY_NAME) %>% 
  mutate(ID = case_when(
    COHORT == "PARKFIT"~ tolower(ID), 
    TRUE~ID
  ))

PhenoFile = inner_join(temp1, Pheno, by=c("COHORT", "PID"="ID"))
write.csv(PhenoFile, "PhenoFile.csv", row.names = F)