rm(list=ls())
library(data.table)
library(dplyr)
imputed = fread("../data/imputedIDs.csv")
imputed %>% with(table(DATASET))
csg = fread("../data/CSGall.csv")
ltg = fread("../data/LTGall.csv")
csg %>% distinct(STUDY_NAME, ID) %>% nrow
ltg %>% distinct(STUDY_NAME, ID) %>% nrow
csg %>% filter(STUDY_NAME !="Oslo") %>% distinct(STUDY_NAME, ID) %>% nrow
# Oslo was not in LSG file (Because baseline visit of some subjects were not able to be estimated)

unique(imputed$DATASET)
# First identify HBS, PDBP, PPMI
d1 = csg %>% filter(STUDY_NAME %in% c("HBS", "PDBP", "PPMI", "Oslo")) %>% select(STUDY_NAME, ID) %>% 
  rename(DATASET = STUDY_NAME,
         PHENOID = ID) %>% 
  mutate(IN = 1, 
         DATASET = ifelse(DATASET=="Oslo", "OSLO", DATASET))



imputed %>% filter(DATASET == "HBS") %>% head  
d1 %>% filter(DATASET == "HBS") %>% head  
# I could not find the HBS IID-PHENOID match file for HBS
csg_HBS = d1 %>% filter(DATASET=="HBS")



imputed %>% filter(DATASET == "PDBP") %>% head  
d1 %>% filter(DATASET == "PDBP") %>% head  
# Likely the first 10 characters are the phenotype ID
csg_PDBP = d1 %>%  filter(DATASET == "PDBP")
imputed_PDBP = imputed %>% filter(DATASET == "PDBP") %>% 
  mutate(PHENOID = substring(IID, 1, 10)) %>% 
  left_join(., csg_PDBP, by =c("PHENOID", "DATASET")) %>% 
  mutate(PHENOID = ifelse(is.na(IN), "NotInPheno", PHENOID))
imputed_PDBP %>% with(table(IN, useNA = "always"))
missing1 = imputed_PDBP %>% filter(!is.na(IN)) %>% 
  anti_join(csg_PDBP, .,  by =c("PHENOID", "DATASET"))


# PPMI 
imputed %>% filter(DATASET == "PPMI") %>% head  
d1 %>% filter(DATASET == "PPMI") %>% head  
# Likely the first 4 characters are the phenotype ID
csg_PPMI = d1 %>%  filter(DATASET == "PPMI")
imputed_PPMI = imputed %>% filter(DATASET == "PPMI") %>% 
  mutate(PHENOID = substring(IID, 1, 4)) %>% 
  left_join(., csg_PPMI, by =c("PHENOID", "DATASET")) %>% 
  mutate(PHENOID = ifelse(is.na(IN), "NotInPheno", PHENOID))
imputed_PPMI %>% with(table(IN, useNA = "always"))
missing2 = imputed_PPMI %>% filter(!is.na(IN)) %>% 
  anti_join(csg_PPMI, .,  by =c("PHENOID", "DATASET"))



# Oslo
imputed %>% filter(DATASET == "OSLO") %>% head  
d1 %>%  filter(DATASET=="OSLO") %>% nrow
# I could not find the matching file...
Oslo = fread("../data/Oslo.csv") %>% select(STUDY_NAME, ID) %>% 
  rename(DATASET = STUDY_NAME,
         PHENOID = ID) %>% 
  mutate(IN = 1, 
         DATASET = ifelse(DATASET=="Oslo", "OSLO", DATASET)) %>% distinct(PHENOID, .keep_all = T)
  


imputed2 = imputed %>% filter(DATASET %in% c("HBS","OSLO")) %>% 
  mutate(PHENOID = "Unknown",
         IN = NA) %>% 
  bind_rows(., imputed_PPMI, imputed_PDBP) %>% 
  mutate(IN = ifelse(is.na(IN), 0, 1))

missings = bind_rows(Oslo, csg_HBS, missing1, missing2)

write.csv(imputed2, "../outputs/IDs_GxP_imputed_match.csv", row.names = F)
write.csv(missings, "../outputs/IDs_GxP_Unknown_IID.csv", row.names = F)
