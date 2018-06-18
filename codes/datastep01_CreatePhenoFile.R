# Phenofile
# Keep IDs available in both ltg_all and imputed files.

setwd("..")
# install.packages("openxlsx")
require(openxlsx)
require(dplyr)
require(data.table)
###################################
# Add phenotype ID column for the imputed individuals
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

###############################
# Read phenotype data
##############################
temp1 =fread('data/LTGall.csv') %>% data.frame %>%
  rename(COHORT = STUDY_NAME) %>% 
  mutate(COHORT= case_when(
    COHORT == "Coriell" ~ "CORIELL",
    COHORT == "PreCEPT" ~ "PRECEPT",
    COHORT == "ParkFit" ~ "PARKFIT",
    COHORT == "ParkWest" ~ "PARKWEST",
    COHORT == "Udall" ~ "UDALL",
    TRUE ~ COHORT
  )) %>% select(-EURO)

Oslo = fread("data/Oslo.csv") %>% 
  mutate(COHORT = "OSLO") %>% 
  mutate(YEARfDIAG = BLDfDIAG + TSTART/365.25) %>% 
  select(names(temp1))

ltg_all = rbind(Oslo, temp1) %>% 
  arrange(COHORT, ID, TSTART)
write.csv(ltg_all, "Longitudinal.csv", row.names = F)
  
#############################
# PhenoFile
#############################
temp1 = fread("ID_imputed.csv") 
temp2 = fread("Longitudinal.csv") %>% 
  mutate(ID = case_when(
    COHORT == "PARKFIT"~ tolower(ID), 
    TRUE~ID
  ))
PhenoFile = inner_join(temp1, temp2, by=c("COHORT", "PID"="ID"))

temp3 = anti_join(temp2, temp1, by=c("COHORT", "ID"="PID")) %>% 
  distinct(COHORT, ID) ## Unmatched in Longitudinal data
temp3 %>% group_by(COHORT) %>% 
  mutate(N = n()) %>% distinct(COHORT, N)

write.csv(temp3, "PhenoInGenoMiss.csv", row.names = F)
