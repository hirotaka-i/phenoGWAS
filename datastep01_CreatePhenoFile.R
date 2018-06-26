# Create phenofile with genetics data
# Use phenofile created from raw data
# HBS data has not genetic-phenotypic table so cannot join.

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
    TRUE~ strsplit(ID, "_")[[1]][2])) %>% data.frame %>%
# This file is tricky so need to correct the ID for combining with fam/bin/bed files
  mutate(ID = case_when(
    DATASET %in% c("CORIELL", "DATATOP", "PARKFIT", "PARKWEST", "PICNICS", "PRECEPT", "SCOPA")~ PID,
    TRUE~ ID)) %>% 
  arrange(DATASET, ID)
temp3 = read.xlsx("data/progression_GWAS_imputed_IDs.xlsx", sheet = 2) %>% 
  arrange(DATASET, ID) %>% 
  rename(PID = Clincal_ID)
temp4 = bind_rows(temp2, temp3) %>% 
  arrange(DATASET, ID) %>% 
  rename(COHORT = DATASET)
write.csv(temp4, "outputs/ID_imputed.csv", row.names = F)

#############################
# PhenoFile with genotyping data
# PhenoFile should be in PDcohorts/codes. Before retrieve, make sure the file is latest
#############################
temp1 = fread("outputs/ID_imputed.csv")
Pheno = fread("data/Pheno13_2018-06-19.csv", stringsAsFactors = F) %>% 
  rename(COHORT = STUDY_NAME) %>% 
  mutate(ID = case_when(
    COHORT == "PARKFIT"~ tolower(ID), 
    TRUE~ID
  ))

PhenoFile = inner_join(temp1, Pheno, by=c("COHORT", "PID"="ID"))
write.csv(PhenoFile, "outputs/PhenoFile.csv", row.names = F)
print("Complete")

rm(list=ls())
