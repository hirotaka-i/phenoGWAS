library(data.table)
library(dplyr)
d1 = fread('../data/ID_Matched.csv')
d2 = fread("../data/IDtoImpute.csv")

data = fread("/Users/iwakih2/Downloads/META/CSGall.csv") %>% select(ID, STUDY_NAME) %>% mutate(IN = 1) %>%
  filter(STUDY_NAME %in% c("Coriell", "DATATOP", 'ParkFit', "ParkWest", "PICNICS", "PreCEPT", "Predict", "SCOPA", "DIGPD")) %>% 
  mutate(DATASET = case_when(
    STUDY_NAME == "Coriell" ~ "CORIELL", 
    STUDY_NAME == "ParkFit"~"PARKFIT",
    STUDY_NAME == "ParkWest"~ "PARKWEST",
    STUDY_NAME == "PreCEPT"~ "PRECEPT", 
    STUDY_NAME == "Predict"~ "PREDICTPD",
    TRUE ~ STUDY_NAME
  ))
d5 = full_join(d2, data, by =c("IID"= "ID", "DATASET"))
d5 %>% with(table(DATASET, !is.na(IN),useNA = "always"))
# if number of TRUE == FALSE, that mean non of the mumbers in the cohort matched. 
# PARKWEST, PICNICS, PRECEPT, PREDICTPD, SCOPA => IID and ID match
# CORIELL, DATATOP, PARKFIT => non-match
# DIGPD => partial match?

# DIGPD 
d2 %>% filter(DATASET=="DIGPD") %>% nrow
data %>% filter(STUDY_NAME=="DIGPD") %>% nrow
# all IDs in phenotyped data were in the givein ID_to_Imputed dataset.
indata = data %>% filter(STUDY_NAME=="DIGPD") %>% select(ID) %>% t %>% as.vector


  
# Create PHENOID for those cohorts with IID = ID
IDs_GxP1=d2 %>% 
  filter(DATASET %in% c("PARKWEST", "PICNICS", "PRECEPT", "PREDICTPD", "SCOPA", "DIGPD")) %>% 
  mutate(PHENOID = IID) %>% 
  mutate(PHENOID = ifelse(DATASET=="DIGPD" & !(IID %in% indata), "NotInPheno", PHENOID))
IDs_GxP1 %>% with(table(DATASET))

# Parkfit and Corieel
d6 = d1 %>% filter(Cohort %in% c("Parkfit", "Coriell")) %>% 
  mutate(DATASET = case_when(
    Cohort == "Coriell" ~ "CORIELL", 
    Cohort == "Parkfit"~"PARKFIT",
    Cohort == "ParkWest"~ "PARKWEST",
    Cohort == "PreCEPT"~ "PRECEPT", 
    Cohort == "Predict"~ "PREDICTPD",
    TRUE ~ Cohort
  ))

d7 = d2 %>% filter(DATASET %in% c("PARKFIT", "CORIELL"))
d8 = left_join(d7, d6, by =c("IID"= "GeneID", "DATASET"))
d8 %>% with(table(DATASET, Flag, useNA = "always"))
# all of PARKFIT and part of Coriel is in d1

# D1 didn't match in the first trial
d8 %>% filter(is.na(Flag)) %>% head
d9 = d8 %>% filter(is.na(Flag)) %>% .[,1:7] %>% 
  mutate(ID = substring(IID, 9, nchar(.))) %>% 
  left_join(., d6, by =c("ID"= "GeneID", "DATASET"))
d9 %>% with(table(DATASET, Flag, useNA = "always"))
# They are matched in the next trial (In the csv file which was genotyped previously)
d8 %>% filter(!is.na(Flag)) %>% head
IDs_GxP2 = d8 %>% filter(!is.na(Flag)) %>% .[,c(1:7, 10)] %>% 
  rename(PHENOID = SubjID)
d9 %>% head
IDs_GxP3 = d9 %>% filter(!is.na(Flag)) %>% .[,c(1:7, 11)] %>% 
  rename(PHENOID = SubjID)


# DATATOP
d2 %>% filter(DATASET=="DATATOP") %>% head
data %>% filter(DATASET=="DATATOP") %>% head
datatopID = fread("/Users/iwakih2/Downloads/PhenoGWAS/data/NallsID_n440.txt")
d10 = d2 %>% filter(DATASET=="DATATOP") %>% 
  left_join(., datatopID, by=c("IID"="NallsID")) %>% 
  mutate(PHENOID = as.character(ID), 
         ID = NULL)
d10 %>% with(table(is.na(PHENOID)))
# All were matched. 
anti_join(d10, data, by = c("PHENOID"="ID"))
# All is in the phenotype data file.
IDs_GxP4 = d10


IDs_GxP = bind_rows(list(IDs_GxP1, IDs_GxP2, IDs_GxP3, IDs_GxP4))
anti_join(d2, IDs_GxP, by=c("DATASET", "IID"))
anti_join(IDs_GxP, d2, by=c("DATASET", "IID"))
# Both are unique now.

write.csv(IDs_GxP, "../outputs/IDs_GxP_2Bimputed_match.csv", row.names = F)
