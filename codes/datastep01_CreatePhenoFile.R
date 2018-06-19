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
# Create phenotype data
##############################
VARS = c("STUDY_NAME", "ID", "TSTART", "DOPA", "AGONIST",
         "HYPOSMIA", "DEMENTIA","MOTORFLUX", "DYSKINESIAS", "DEPR", "RL","CONST", "RBD", "SLEEP", "INS", "HY3",
         "HY", "UPDRS1_scaled", "UPDRS2_scaled", "UPDRS3_scaled", "UPDRS4_scaled", "UPDRS1","UPDRS2","UPDRS3", "UPDRS4", 
         "oldUPDRS","MDS_UPDRS", "MMSE", "MOCA", "SEADL",
         "AGEatBL","FEMALE","FAMILY_HISTORY", "YEARSEDUC","BLDfDIAG","AAO","AD")

# CIRIELL
#######################################################
## DEPR
Coriell_Depr = fread('../PDcohorts/Coriell/formatted/bdi.csv') %>% 
  mutate(DEPR = ifelse(BDI_total>14, 1, 0)) %>% 
  select(SUBJID, daysB, DEPR) 
## DEMO
Coriell_Demo = fread('../PDcohorts/Coriell/formatted/demo.csv') %>% 
  rename(AGEatBL = age,
         EUROPEAN = RAWHITE) %>% 
  mutate(FEMALE = 2-gender,
         YEARSEDUC = case_when(EDUC==1 ~ 9,
                               EDUC==2 | EDUC==3 ~ 12,
                               EDUC==4 | EDUC==5 ~ 14, 
                               EDUC==6 ~ 16, 
                               EDUC==7 ~ 18)) %>% 
  select(SUBJID, FEMALE, EUROPEAN, YEARSEDUC, AGEatBL) 

Coriell_Death = fread('../PDcohorts/Coriell/formatted/status.csv') %>% 
  select(SUBJID, death, days_died) 

## DIAGFEAT only at Baseline and the study ends
Coriell_Diag = fread('../PDcohorts/Coriell/formatted/diagfeat.csv')

## Family History
Coriell_FH = fread('../PDcohorts/Coriell/formatted/famhxpd.csv') %>% 
  data.frame %>% 
  mutate(FAMILY_HISTORY = rowSums(.[,c("biodadpd", "biomompd", "fulsibpd")], na.rm = T)) %>% 
  mutate(FAMILY_HISTORY = ifelse(FAMILY_HISTORY>0, 1, 0)) %>% 
  distinct(SUBJID, .keep_all=T) %>% 
  select(SUBJID, FAMILY_HISTORY)

## anti-Parkisonism medication
Coriell_Drug = fread('../PDcohorts/Coriell/formatted/led_daily.csv') %>% 
  mutate(DOPA = if_else(Levodopa>0, 1, 0),
         AGONIST = if_else(Dopamine>0, 1, 0))%>% 
  rename(daysB = daysb) %>% 
  select(SUBJID, daysB, DOPA, AGONIST)

## SEADL
Coriell_SEADL = fread('../PDcohorts/Coriell/formatted/modseadl.csv') %>% 
  rename(SEADL = mseadlg) %>% 
  select(SUBJID, daysB, SEADL)

## Diseaes duration need to calculate for AAO later
Coriell_AAO = fread('../PDcohorts/Coriell/formatted/pdfeat.csv') %>% 
  mutate(BLDfDIAG = duration - daysB/365.25) %>% 
  select(SUBJID, BLDfDIAG)

## PDQ39
Coriell_PDQ39 = fread('../PDcohorts/Coriell/formatted/pdq39.csv') %>% 
  rename(PDQ39 = PDQ39_ADL) %>% 
  select(SUBJID, daysB, PDQ39)

## SCOPACOG
Coriell_Dementia = fread('../PDcohorts/Coriell/formatted/scopacog.csv') %>% 
  mutate(DEMENTIA = ifelse(scopacog<23, 1, 0)) %>% 
  select(SUBJID, daysB, DEMENTIA)

## UPDRS
Coriell_UPDRS123 = fread('../PDcohorts/Coriell/formatted/updrs.csv') %>% 
  rename(UPDRS1 = UPDRS_mental,
         UPDRS2 = UPDRS_ADL, 
         UPDRS3 = UPDRS_motor,
         oldUPDRS = UPDRS_total) %>%
  mutate(HY3 = ifelse(p3pstabl>1, 1, 0)) %>% 
  select(SUBJID, daysB, paste('UPDRS', 1:3, sep = ''), oldUPDRS, HY3)

Coriell_UPDRS4 = fread('../PDcohorts/Coriell/formatted/updrs4.csv') %>% 
  data.frame %>% 
  mutate(UPDRS4 = rowSums(.[, grep("^p4", names(.))], na.rm=T),
         DYSKINESIAS = ifelse(p4adurat>0, 1, 0)) %>% 
  mutate(MOTORFLUX = ifelse(p4boff>0 | p4adurat>0, 1, 0)) %>% 
  select(SUBJID, daysB, UPDRS4, MOTORFLUX, DYSKINESIAS)

Coriell_UPDRS = left_join(Coriell_UPDRS123, Coriell_UPDRS4, by = c('SUBJID', 'daysB')) %>% 
  mutate(oldUPDRS = oldUPDRS + UPDRS4)

base = Coriell_UPDRS %>% 
  arrange(daysB) %>% 
  distinct(SUBJID, .keep_all=T) %>% 
  summarize_at(vars(contains('UPDRS')), funs(mean(., na.rm=T), sd(., na.rm=T)))

Coriell_UPDRS = Coriell_UPDRS %>% 
  mutate(UPDRS1_scaled = (UPDRS1 - base$UPDRS1_mean)/base$UPDRS1_sd,
         UPDRS2_scaled = (UPDRS2 - base$UPDRS2_mean)/base$UPDRS2_sd,
         UPDRS3_scaled = (UPDRS3 - base$UPDRS3_mean)/base$UPDRS3_sd,
         UPDRS4_scaled = as.numeric(scale(UPDRS4)),
         UPDRS_scaled = (oldUPDRS - base$oldUPDRS_mean)/base$oldUPDRS_sd)

# Make a longitudinal database
Coriell = full_join(Coriell_Demo, Coriell_FH, by = 'SUBJID') %>% 
  full_join(., Coriell_AAO, by = 'SUBJID') %>% 
  left_join(., Coriell_Death, by = 'SUBJID') %>% 
  inner_join(Coriell_UPDRS, ., by = 'SUBJID') %>% 
  full_join(., Coriell_Dementia, by=c('SUBJID', 'daysB')) %>%
  full_join(., Coriell_Depr, by=c('SUBJID', 'daysB')) %>% 
  full_join(., Coriell_PDQ39, by=c('SUBJID', 'daysB')) %>% 
  full_join(., Coriell_SEADL, by=c('SUBJID', 'daysB')) %>% 
  left_join(., Coriell_Drug, by=c('SUBJID', 'daysB'))

Coriell = Coriell %>% 
  mutate(STUDY_NAME = 'CORIELL', 
         AAO = AGEatBL - BLDfDIAG, MDS_UPDRS=NA,
         AD = AGEatBL + days_died/365.25, death =NULL, days_died=NULL,
         HYPOSMIA = NA, RL = NA, HT = NA, CONST = NA, RBD = NA, SLEEP = NA, INS = NA,
         HY = NA, MMSE = NA, MOCA=NA) %>% 
  rename(TSTART = daysB)



# identify the ones with TSTART <0 and give locf value for more than 2 screening.
Check_ID = Coriell %>% 
  filter(TSTART<0) %>% 
  distinct(SUBJID)

Check = inner_join(Check_ID, Coriell, by='SUBJID') %>% 
  filter(TSTART <= 0) %>% 
  arrange(TSTART) %>% 
  group_by(SUBJID) %>% 
  mutate_all(funs(na.locf(., na.rm=F))) %>% 
  arrange(desc(TSTART)) %>% 
  distinct(SUBJID, .keep_all=T) %>% 
  data.frame

# Merge and reset TSTART from 0.
CORIELL = Coriell %>% 
  filter(TSTART >= 0) %>% 
  rbind(., Check) %>% 
  arrange(SUBJID, TSTART) %>% 
  group_by(SUBJID) %>% 
  mutate(TSTART = TSTART - first(TSTART)) %>% 
  data.frame %>% 
  rename(ID = SUBJID) %>% 
  select(VARS)


# DATATOP 
############################################################
lt = fread("../PDcohorts/DATATOP/DATATOP_longitudinal.txt", header = T) %>% 
  rename(DATE = VISIT_DAYS,
         oldUPDRS = UPDRS, 
         HYPOSMIA = ANOSMIA,
         HYPOSMIA_ONSETDAYS = ANOSMIA_ONSETDAYS) %>% 
  arrange(ID, DATE) %>% 
  distinct(ID, DATE, .keep_all = T) %>% 
  data.frame

lt_TMP = lt[, c(3:8, 21:23)]
names(lt_TMP) # Only continous
# Binomial should be changed into an ordinary long format after adding the event date observation
TMP_VARS= c("HYPOSMIA", "DYSKINESIAS", "DEPR", "CONST", "INS", "SLEEP") # Binom Vars
for (i in 1:length(TMP_VARS)){
  lt_TMP2 = lt[, c(3,4, grep(TMP_VARS[i], names(lt)))]
  print("Before")
  lt_TMP2 %>% filter(!is.na(lt_TMP2[,4])) %>% head(n=20) %>% print
  lt_TMP3 = lt_TMP2 %>% filter(!is.na(lt_TMP2[,4])) %>%  
    distinct(ID, .keep_all = T) %>% 
    mutate(DATE = ifelse(.[,2] <= .[,4], .[,4], .[,2])) %>% 
    bind_rows(lt_TMP2, .) %>% 
    arrange(ID, DATE) %>% 
    distinct(ID, DATE, .keep_all = T)
  
  lt_TMP3[,3] = ifelse(!is.na(lt_TMP3[,4]) & lt_TMP3$DATE < lt_TMP3[,4], 0, lt_TMP3[,3])
  print("After")
  lt_TMP3 %>% filter(!is.na(.[,4])) %>% head(n=20) %>% print
  # GOOD. Marge this file with the original file
  lt_TMP = full_join(lt_TMP, lt_TMP3[,1:3], by= c("ID", "DATE")) %>% 
    arrange(ID, DATE) %>%
    distinct(ID, DATE, .keep_all = T)
}

# Give LED, DOPA, and AGONIT values for new event date using locf value
lt_TMP = lt_TMP %>% 
  mutate_at(vars("DOPA", "AGONIST", "LED"), funs(na.locf(., na.rm = F)))


lt = lt_TMP %>% 
  mutate(DEMENTIA= if_else(MMSE<27, 1, 0),
         MOTORFLUX = NA,
         RL = NA,
         HT = NA,
         RBD = NA,
         UPDRS1 = NA, 
         UPDRS2 = NA, 
         UPDRS3 = NA,
         UPDRS4 = NA,
         MDS_UPDRS = NA, 
         MOCA = NA,
         PDQ39 = NA, 
         HY3 = ifelse(HY>=3.0, 1, 0)) %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE)))%>%  # DAys from the first visit
  data.frame
##################################

###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (oldUPDRS - mean(Baseline$oldUPDRS, na.rm=T))/sd(Baseline$oldUPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4)))


cs =  fread("../PDcohorts/DATATOP/DATATOP_crosssectional.txt", header = T) %>%
  mutate(AGEatBL = AGE,
         BLDfDIAG = AGE-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY== 1, 1, 0)) %>% 
  distinct(ID, .keep_all=T)

DATATOP = cs %>% 
  inner_join(lt, ., by = 'ID') %>% 
  select(VARS)


#DIGPD
######################################################################################################
lt = fread("../PDcohorts/DIGPD/DIGPD_longitudinal.txt", header = T) %>% 
  mutate(DATE=as.Date(VISIT, format='%m/%d/%Y')) %>% 
  rename(UPDRS1 = MDS_UPDRS_Subscale_1, 
         UPDRS2 = MDS_UPDRS_Subscale_2, 
         UPDRS3 = MDS_UPDRS_Subscale_3, 
         UPDRS4 = MDS_UPDRS_Subscale_4) %>% 
  mutate(MDS_UPDRS = rowSums(.[,grep('UPDRS*', colnames(.))], na.rm = T),
         oldUPDRS = NA,
         RBD=NA, 
         ID = UniqID, 
         DYSKINESIAS = ifelse(DYSKINESIAS ==2, NA, DYSKINESIAS),
         DOPA = ifelse(DOPA==11, NA, DOPA),
         HY3 = ifelse(HY>=3.0, 1, 0)) %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame


# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (1 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

lt = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() == 1) %>% 
  mutate(idx=1) %>% data.frame %>% #idx == number of the observation in the same ID, Date
  rbind(., lt_process_duplicate)


###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (MDS_UPDRS - mean(Baseline$MDS_UPDRS, na.rm=T))/sd(Baseline$MDS_UPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4)))

cs =  fread("../PDcohorts/DIGPD/DIGPD_crosssectional.txt", header = T) %>%
  mutate(ID = UniqID,
         AGEatBL = AGE,
         BLDfDIAG = AGE-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY=='European', 1, 0))

DIGPD = inner_join(lt, cs, by = 'ID') %>% mutate(STUDY_NAME = "DIGPD") %>% select(VARS)


# HBS
#####################################################################
lt = fread("../PDcohorts/HBS/HBS_longitudinal.csv", header = T) %>% 
  rename(UPDRS1 = UPDRS_Subscale_1, 
         UPDRS2 = UPDRS_Subscale_2,
         UPDRS3 = UPDRS_Subscale_3,
         UPDRS4 = UPDRS_Subscale_4) %>% 
  mutate(oldUPDRS = UPDRS1 + UPDRS2 + UPDRS3 + UPDRS4,
         DATE = as.Date(VISIT, '%m/%d/%Y'), 
         RBD = NA, 
         DEMENTIA = if_else(MMSE < 27, 1, 0),
         HY3 = ifelse(HY>=3, 1, 0)) %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  distinct(ID, TSTART, .keep_all=T) %>% 
  data.frame 

###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (oldUPDRS - mean(Baseline$oldUPDRS, na.rm=T))/sd(Baseline$oldUPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4))) 


cs =  fread("../PDcohorts/HBS/HBS_crosssectional.tab", header = T) %>% 
  mutate(AGEatBL = AGE, 
         BLDfDIAG = AGEatBL - AAO, # Years from diagnosis to baseline
         EUROPEAN = if_else(ANCESTRY=='European', 1, 0)) %>% 
  distinct(ID, .keep_all=T)

HBS = left_join(lt, cs, by = 'ID') %>% mutate(STUDY_NAME = "HBS") %>%  select(VARS)


# OSLO
# The data has both retrospective and prospective data
# HY3, HY, UPDRS. No baseline outcomes.
##########################################################
#####

d1 = fread("../PDcohorts/Oslo/Oslo_UPDRS3_scorings.txt", header = T) %>% 
  mutate(DAYSfDIAG = days_from_Jan1st_YoD,
         UPDRSnum=1) %>% 
  filter(!is.na(DAYSfDIAG)) %>% 
  filter(!is.na(UPDRS3_ON)) %>% 
  arrange(DAYSfDIAG) %>% 
  group_by(ID) %>% 
  mutate(UPDRSnum = cumsum(UPDRSnum)) %>% 
  data.frame %>% arrange(ID) %>% 
  select(-days_from_Jan1st_YoD) %>% 
  rename(UPDRS3 = UPDRS3_ON)

Baseline = d1 %>% filter(UPDRSnum==1) 
UPDRS =  d1 %>% 
  mutate(UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T))

rm(d1, Baseline)

HY = fread("../PDcohorts/Oslo/Oslo_HoehnYahr_scorings.txt", header = T) %>% 
  mutate(DAYSfDIAG = days_from_Jan1st_YoD,
         HYnum=1) %>% 
  filter(!is.na(DAYSfDIAG)) %>% 
  filter(!is.na(Hoehn_Yahr_score)) %>% 
  arrange(DAYSfDIAG) %>% 
  group_by(ID) %>% 
  mutate(HYnum = cumsum(HYnum)) %>% 
  data.frame %>% arrange(ID) %>% 
  select(-days_from_Jan1st_YoD) %>% 
  rename(HY = Hoehn_Yahr_score) # Do not define HY3 here. Use other provided data for HY3

HYUPDRS = full_join(UPDRS, HY, by=c("ID", "DAYSfDIAG")) %>% arrange(ID, DAYSfDIAG)
HYUPDRS %>%  group_by(ID, DAYSfDIAG) %>% mutate(N = n()) %>% filter(N>1)
HYUPDRS %>% filter(ID == "PDRH049")
HYUPDRS %>% filter(ID == "PDRH175")
HYUPDRS %>% filter(ID == "PDRH242")

# Take the larger value if observed on the same date.
HYUPDRS2 = HYUPDRS %>% arrange(ID, DAYSfDIAG, desc(UPDRS3), desc(HY)) %>% distinct(ID, DAYSfDIAG, .keep_all = T)

d2 = fread("../PDcohorts/Oslo/Oslo_genotypes_data_PCs.txt", header = T, na.strings = "-9") %>% 
  mutate(ID = IID, 
         FEMALE = SEX-1,
         EUROPEAN = ethnicity,
         YEARSEDUC = NA,
         AGEatBL = age_at_study, 
         BLDfDIAG = age_at_study - age_at_diagnosis,
         AAO = age_at_diagnosis, 
         AD = ifelse(dead==1, days_from_diagnosis_dead_or_censored/365.25 + age_at_diagnosis, NA), 
         FAMILY_HISTORY = ifelse(family_history_1st_degree_relative==1, 1, 0))

HY3_Base = d2 %>% filter(!is.na(progressed_to_HY3)) %>% 
  mutate(DAYSfDIAG=0, HY3 = 0) %>% 
  select(ID, DAYSfDIAG, HY3)
HY3_End  = d2 %>% filter(!is.na(progressed_to_HY3)) %>% 
  mutate(DAYSfDIAG = days_from_diagnosis_HY3_or_censored, HY3 = progressed_to_HY3) %>%
  select(ID, DAYSfDIAG, HY3)
HY3 = bind_rows(HY3_Base, HY3_End)

# HYUPDRS contains retrospective information (HY and UPDRS before the prospective cohort enrollment)
temp = left_join(d2, HYUPDRS2, by = "ID") %>% arrange(DAYSfDIAG) %>%
  mutate(YEARfDIAG = DAYSfDIAG/365.25) %>% distinct(ID, .keep_all = T) %>% 
  select(ID, AAO, AGEatBL,BLDfDIAG, YEARfDIAG) %>% 
  mutate(FIRST_FOLLOW_UP = YEARfDIAG-BLDfDIAG) %>% arrange(FIRST_FOLLOW_UP) %>% filter(FIRST_FOLLOW_UP<0)
print(paste(nrow(temp), "observations are started before the entry (retrospecitve)"))
temp %>% head


# People in d3 but not in HYUPDRS2 -> 5 observation
anti_join(d2, HYUPDRS2, by = "ID") %>% select(ID)
# People in HYUPDRS2 but not in d3 = 188
anti_join(HYUPDRS2, d2, by = "ID") %>% select(ID) %>% nrow

# Longitudinal data
OSLO_LT = full_join(HY3, HYUPDRS2, by = c("ID", "DAYSfDIAG")) %>% 
  left_join(d2, ., by="ID") %>% 
  arrange(DAYSfDIAG) %>% 
  mutate(TSTART = DAYSfDIAG - BLDfDIAG*365.25) %>% # TSTART is negative for retrospective data
  mutate(STUDY_NAME="OSLO",
         DOPA=NA, AGONIST=NA, HYPOSMIA=NA, DEMENTIA=NA, MOTORFLUX=NA, DYSKINESIAS=NA, 
         DEPR=NA, RL=NA, HT=NA, CONST=NA, RBD = NA, SLEEP=NA, INS=NA, UPDRS1_scaled = NA,
         UPDRS2_scaled=NA, UPDRS4_scaled=NA, UPDRS_scaled =NA, MMSE =NA, MOCA=NA, SEADL=NA, 
         PDQ39=NA, UPDRS1=NA, UPDRS2=NA, UPDRS4=NA, oldUPDRS=NA, MDS_UPDRS=NA) %>% 
  arrange(ID, TSTART) %>% 
  select(VARS)
anti_join(OSLO1, OSLO, by = "ID") %>% select(ID)

OSLO_LT %>% with(hist(TSTART))
OSLO_LT %>% filter(-180<TSTART & TSTART<180) %>% distinct(ID) %>% nrow


# Parkfit
###################################################################################
Parkfit_base = fread('../PDcohorts/Parkfit/Parkfit_base.csv', na.strings = "999") %>%
  mutate(TSTART = 0,
         DEMENTIA = ifelse(MMSE<27, 1, 0)) %>% 
  rename(SUBJID = PTNR,
         AGEatBL = Age,
         HY = HY_0,
         BLDfDIAG = Duration_yrs,
         FEMALE = Sex,
         oldUPDRS = UPDRS_0) %>% data.frame

vars = c('SUBJID', 'TSTART', 'HY', 'oldUPDRS', 'MMSE', 'DEMENTIA')
Parkfit_base_f = Parkfit_base %>% select(vars)
Parkfit_basic = Parkfit_base %>% select(SUBJID, FEMALE, AGEatBL, BLDfDIAG)

Parkfit_res  = fread('../PDcohorts/Parkfit/24mnths_PFdata_clinical characteristics.csv') %>% 
  mutate_if(is.factor, as.character) %>% 
  rename(SUBJID = Patientnumber) %>%
  mutate_at(vars(names(.)[2:9]), as.numeric) %>% 
  na_if(., 999) %>% 
  mutate(TSTART = 720, 
         MMSE = NA,
         DEMENTIA = NA) %>% 
  rename(HY = HY_24mths,
         oldUPDRS = UPDRS_24mnths) 

PARKFIT = Parkfit_res%>% 
  select(vars) %>% 
  rbind(., Parkfit_base_f) %>% 
  left_join(., Parkfit_basic, by='SUBJID') %>% 
  rename(ID = SUBJID) %>% 
  mutate(AAO = AGEatBL - BLDfDIAG, 
         AD = NA, FAMILY_HISTORY=NA, 
         EUROPEAN=NA, YEARSEDUC=NA, DOPA=NA, AGONIST=NA,
         HYPOSMIA = NA, MOTORFLUX=NA, DYSKINESIAS=NA, DEPR=NA, RL=NA,
         HY3=ifelse(HY>=3, 1, 0),
         UPDRS1_scaled=NA, UPDRS2_scaled=NA, UPDRS3_scaled=NA, UPDRS4_scaled=NA,
         UPDRS_scaled = (oldUPDRS - mean(Parkfit_base_f$oldUPDRS, na.rm=T))/sd(Parkfit_base_f$oldUPDRS, na.rm = T),
         UPDRS1=NA, UPDRS2=NA, UPDRS3=NA, UPDRS4=NA, MDS_UPDRS=NA, 
         HT=NA, CONST=NA, RBD=NA, SLEEP=NA, INS=NA, 
         MOCA=NA, SEADL=NA, PDQ39=NA, 
         STUDY_NAME = 'PARKFIT') %>% 
  select(VARS)


#######################################################
data <- read.table("../PDcohorts/ParkWest/working.delim", header = T, sep = ";", stringsAsFactors = F)
dat <- subset(data, V9_Visit != " " & ID != " " & Study_name != " ")
dat$V9_CONSTIPAT <- ifelse(is.na(dat$V9_CONSTIPAT_OFF) & is.na(dat$V9_CONSTIPAT_ON), NA,
                           ifelse((!is.na(dat$V9_CONSTIPAT_OFF) & dat$V9_CONSTIPAT_OFF == 1) | 
                                    !is.na(dat$V9_CONSTIPAT_ON) & dat$V9_CONSTIPAT_ON == 1, 1, 0))

BL <- dat[,c("Study_name","STUDY_SITE","ID","FEMALE","CASE","FAMILY_HISTORY","Age_Screen", "Age_BL", "Age_DEATH","YEARSEDUC","BL_Visit","BL_UPO_part1_tot","BL_UPO_part2_tot","BL_UPO_part3_tot","BL_UPO_part4_tot","BL_HY","BL_MMSE_tot","BL_GGHC","BL_TGHC","BL_SEADL","BL_HYPOSMIA","BL_DYSKINESIAS","BL_DEPR","BL_MOTORFLUX","BL_CONSTIPAT","BL_INS","BL_SLEEP")]
V9 <- dat[,c("Study_name","STUDY_SITE","ID","FEMALE","CASE","FAMILY_HISTORY","Age_Screen", "Age_BL", "Age_DEATH","YEARSEDUC","V9_Visit","V9_UPO_part1_tot","V9_UPO_part2_ON_tot","V9_UPO_part3_tot","V9_UPO_part4_tot","V9_HY","V9_MMSE_tot","V9_GGHC","V9_TGHC","V9_SEADL","V9_HYPOSMIA","V9_DYSKINESIAS","V9_DEPR","V9_MOTORFLUX","V9_CONSTIPAT","V9_INS","V9_SLEEP", "V9_agonist", "V9_DPOA")]
names(BL) <- c("STUDY_NAME","STUDY_SITE","ID","FEMALE","CASE","FAMILY_HISTORY","AAO","AGEatBL","AD","YEARSEDUC","VISIT","UPDRS_Subscale_1","UPDRS_Subscale_2","UPDRS_Subscale_3","UPDRS_Subscale_4","HY","MMSE","GGHC","TGHC","SEADL","HYPOSMIA","DYSKINESIAS","DEPR","MOTORFLUX","CONSTIPAT","INS","SLEEP")
BL$DOPA = 0
BL$AGONIST = 0
names(V9) <- c("STUDY_NAME","STUDY_SITE","ID","FEMALE","CASE","FAMILY_HISTORY","AAO", "AGEatBL","AD","YEARSEDUC","VISIT","UPDRS_Subscale_1","UPDRS_Subscale_2","UPDRS_Subscale_3","UPDRS_Subscale_4","HY","MMSE","GGHC","TGHC","SEADL","HYPOSMIA","DYSKINESIAS","DEPR","MOTORFLUX","CONSTIPAT","INS","SLEEP", "DOPA", "AGONIST")

combo <- rbind(BL,V9)

lt = combo %>% 
  mutate(DATE=as.Date(VISIT, format="%m/%d/%Y")) %>% 
  rename(UPDRS1 = UPDRS_Subscale_1, 
         UPDRS2 = UPDRS_Subscale_2, 
         UPDRS3 = UPDRS_Subscale_3,
         UPDRS4 = UPDRS_Subscale_4,
         CONST = CONSTIPAT) %>% 
  mutate(STUDY_NAME="PARKWEST", 
         oldUPDRS = rowSums(.[,grep('UPDRS*', colnames(.))], na.rm=T),
         MDS_UPDRS = NA,
         DEMENTIA = if_else(MMSE<27, 1, 0),
         RL = NA,
         HT = NA, 
         RBD = NA, 
         MOCA = NA,
         PDQ39 = NA,
         SLEEP = if_else(SLEEP>9, 1, 0),
         HY3 = ifelse(HY>=3.0, 1, 0)) %>% 
  mutate(BLDfDIAG = AGEatBL - AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = 1, # All european
         FAMILY_HISTORY = ifelse(is.na(FAMILY_HISTORY), 0, ifelse(FAMILY_HISTORY==1, 1, 0))) %>% # 1 missing FAMILY_HISTORY->replace with 0
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame

# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (0 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (oldUPDRS - mean(Baseline$oldUPDRS, na.rm=T))/sd(Baseline$oldUPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4)))

PARKWEST = lt %>% select(VARS)



# PDBP
##############################################################
updrs_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_MDS_UPDRS_2018-02-12T15-37-26.csv", header = T)
updrs1 <- updrs_temp[, c(6, 5, 4, 22, 18, 19, 25, 39, 80, 79, 96, 81, 89, 9)]
colnames(updrs1) <- c("ID", "VISIT", "EVENT",  "CONST", "INS", "SLEEP1", "UPDRS1", "UPDRS2",
                      "UPDRS3", "HY", "UPDRS4", "DYSKINESIAS", "MOTORFLUX", "Age_M")
updrs1$MOTORFLUX <- ifelse(updrs1$MOTORFLUX > 0, 1, 0)
updrs1$DYSKINESIAS =  ifelse(updrs1$DYSKINESIAS > 0, 1, 0)
updrs1$EVENT <- as.character(updrs1$EVENT)

updrs1$CONST <- ifelse(updrs1$CONST > 0, 1, 0)
updrs1$INS <- ifelse(updrs1$INS > 0, 1, 0)
updrs1$SLEEP1 <- ifelse(updrs1$SLEEP1 > 0, 1, 0)

usuful_func = function(data){
  # Return One obs per Visit Type. 
  # Give NA for seemingly wrong VisitDate
  data = data %>% 
    mutate(ID = as.character(ID),
           EVENT = as.character(EVENT),
           VISIT = as.Date(VISIT, format ="%Y-%m-%d")) %>% 
    mutate(MON = ifelse(EVENT=="Baseline", "0 month", EVENT)) %>% 
    mutate(MON = as.numeric(substr(MON, 1, 2)))
  data_dup = data %>% group_by(ID, MON) %>% 
    mutate(dup=n()) %>% data.frame
  print("Duplication per Visit Type (dup>1)")
  data_dup %>% 
    with(table(MON, dup)) %>% print
  # the latest one will be kept with replacement of data from the old one if NA.
  ans = data_dup %>% filter(dup>1) %>% 
    arrange(VISIT) %>% 
    group_by(ID, MON) %>% 
    mutate_all(funs(na.locf(.,na.rm = FALSE))) %>% 
    data.frame %>% 
    arrange(desc(VISIT)) %>% 
    distinct(ID, .keep_all = T) %>% 
    rbind(., subset(data_dup, dup==1)) %>% 
    mutate(dup=NULL) %>% arrange(ID, VISIT)
  # check the inagreement of EVENT order(MON) and Visit date order
  ans_check = ans %>% arrange(MON) %>% 
    group_by(ID) %>% 
    mutate(diff = VISIT - lag(VISIT, default = first(VISIT))) %>% data.frame
  print("Seemingly wrong either VisitDate or VisitTyp in data")
  ans_check %>% filter(diff<0) %>% select(ID, VISIT, EVENT) %>% 
    print
  return = ans_check %>% 
    mutate(VISIT = ifelse(diff>=0, VISIT, NA),
           diff = NULL)
  return(return)
}
updrs1 = usuful_func(updrs1)


### removing trouble makers:
#PDAE054HMM : hoen and yar stage 5 but 0 everything else, all NAs for second entry -> remove
plot(updrs1$UPDRS3, updrs1$HY)
updrs1 <- updrs1[updrs1$ID != "PDAE054HMM",]



#### MOCA and dementia dataset
moca_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_MoCA_2018-02-12T15-35-49.csv", header = T)
moca <- moca_temp[, c(6, 5, 4, 24, 9)]
colnames(moca) <- c("ID", "VISIT", "EVENT",  "MOCA", "Age_M")
moca$DEMENTIA <- ifelse(moca$MOCA < 24, 1, 0)

moca = usuful_func(moca)


##### SEADL
seadl_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_ModSchwabAndEnglandScale_2018-02-12T15-36-55.csv", header = T)
seadl <- seadl_temp[, c(6, 5, 4, 10, 9)]
colnames(seadl) <- c("ID", "VISIT", "EVENT",  "SEADL", "Age_M")
seadl = usuful_func(seadl) %>% 
  mutate(SEADL = SEADL*100)


##### HYPOSMIA
#### no categorization of hyposimia SDORE <= 21 (/40) is defined as HYPOSMIA
upsit_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_UnivOfPennSmellIdenTest_2018-02-12T15-35-36.csv", header = T)
upsit <- upsit_temp[, c(6, 5, 4, 15, 9)]
colnames(upsit) <- c("ID", "VISIT", "EVENT",  "SCORE", "Age_M")
upsit$HYPOSMIA = ifelse(upsit$SCORE <= 21, 1, 0)
upsit$SCORE = NULL
upsit = usuful_func(upsit)



# DEPR
### no GDS for regular data, only for legacy data
dep_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_HDRS_2018-02-12T15-36-42.csv", header = T)
dep_temp <- dep_temp[!is.na(dep_temp$HDRS.HDRS.HDRSTotScore),]
dep <- dep_temp[, c(6, 5, 4, 28, 9)]
colnames(dep) <- c("ID", "VISIT", "EVENT",  "DEPR", "Age_M")
dep$DEPR <- ifelse(dep$DEPR >9,1,0) # Mov Disord. 2007 June 15; 22(8): 1077–1092
dep = usuful_func(dep)

# RLS RBDSQ: RLS(==1), RBD (RBDSQ>5)
sleep_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_MayoSleepQuestionnaire_2018-02-12T15-37-09.csv", header = T)
names(sleep_temp)=c("a", "b", "c")
sleep <- sleep_temp[, c(6, 5, 4, 22, 14, 9)]

colnames(sleep) <- c("ID", "VISIT", "EVENT",  "RL", "RBD", "Age_M")

sleep$RL <- ifelse(sleep$RL == 'Yes', 1, 0)
sleep$RBD <- ifelse(sleep$RBD == "Yes", 1, 0)

sleep = usuful_func(sleep)

# Daytime sleepiness (EPworth, score>=10)
ep_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_EpworthSleepinessScale_2018-02-12T15-36-28.csv", header = T)
ep <- ep_temp[, c(6, 5, 4, 18, 9)]
colnames(ep) <- c("ID", "VISIT", "EVENT",  "SLEEP", "Age_M")

ep$SLEEP <- ifelse(ep$SLEEP > 9, 1, 0)

ep = usuful_func(ep)


# Levodopa or Dopamin Agonist
med_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_PriorAndConcomitantMeds_2018-02-12T15-36-09.csv", header = T)
med_temp <- med_temp[!is.na(med_temp$Study.ID),]
med <- med_temp[, c(6, 5, 4, 10, 9)]
colnames(med) <- c("ID", "VISIT", "EVENT",  "MED", "Age_M")
med$MED <- as.factor(med$MED)
med <- med[med$MED !="",]

med$DOPA <- ifelse(med$MED =="Carbidopa/levodopa orally disintegrating (Parcopa)" |
                     med$MED == "Carbidopa/levodopa/entacapone (Stalevo)" |
                     med$MED == "Carbidopa/levodopa (Atamet, Sinemet)" | 
                     med$MED == "Carbidopa/levodopa sustained release (Sinemet CR)"|
                     med$MED == "Entacapone (Comtan)",1,0)


med$AGONIST <- ifelse(med$MED =="Apomorphine (Apokyn)" | 
                        med$MED == "Bromocriptine (Parlodel)" |
                        med$MED == "Pramipexole (Mirapex)" |
                        med$MED =="Pramipexole extended release (Mirapex ER)" | 
                        med$MED =="Ropinirole (Requip)" | 
                        med$MED =="Ropinirole extended release (Requip XL)" |
                        med$MED =="Rotigotine transdermal patch (Neupro)" ,1,0  )

med <- med[, c(1,2,3,5,6, 7)]

med = usuful_func(med)

# PDQ39
pdq39_temp = fread("../PDcohorts/PDBP/query_result_PDQ39_2016-10-26T15-04-27.csv", header = T)
names(pdq39_temp)
pdq39 = pdq39_temp[,c(6,5,4,50)]
names(pdq39)=c("ID", "VISIT", "EVENT","PDQ39")
pdq39 = pdq39 %>% 
  mutate(PDQ39 = as.numeric(PDQ39)) %>% #ADL
  select(ID, VISIT, EVENT, PDQ39) %>% 
  usuful_func(.)

pdq39 %>% with(hist(PDQ39))

#### combine
combine = 
  full_join(updrs1, dep, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., ep, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., med, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., moca, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., seadl, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., sleep, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., upsit, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(Age_M = ifelse(is.na(Age_M.x), Age_M.y, Age_M.x), 
         Age_M.x=NULL, Age_M.y=NULL) %>% 
  full_join(., pdq39, by = c("ID","EVENT", "MON")) %>% 
  mutate(VISIT = ifelse(is.na(VISIT.x), VISIT.y, VISIT.x), 
         VISIT.x=NULL, VISIT.y=NULL) %>% 
  mutate(VISIT = as.Date(VISIT),
         # MEDICATION for visits
         DOPA = ifelse(is.na(DOPA) & MON==0, 0, DOPA),
         AGONIST = ifelse(is.na(AGONIST) & MON==0, 0, AGONIST)) %>% 
  arrange(MON) %>% 
  group_by(ID) %>% 
  mutate(DOPA = na.locf(DOPA, na.rm=F),
         AGONIST = na.locf(AGONIST, na.rm = F)) %>% 
  data.frame


combine = usuful_func(combine)
# 17 Un defined Visit date -> Put +6 month from the previous visits
combine %>% filter(is.na(VISIT)) %>% 
  select(ID, MON, VISIT, EVENT)

IDs = combine %>% filter(is.na(VISIT)) %>%
  select(ID)

combine1 = combine %>% 
  inner_join(., IDs, by = "ID") %>% 
  arrange(MON) %>% 
  group_by(ID) %>% 
  mutate(VISIT= ifelse(is.na(VISIT), lag(VISIT + 180), VISIT)) %>% 
  data.frame %>% 
  rbind(., combine) %>% 
  distinct(ID, MON, .keep_all = T) %>% 
  arrange(ID, MON)

# With this method, Visit is OK but Age doesn't agree with visits.
combine1 %>% arrange(MON) %>% 
  group_by(ID) %>% 
  mutate(diff = Age_M - lag(Age_M, default = first(Age_M))) %>%
  filter(diff<0) %>% select(ID, VISIT, EVENT, Age_M)


# Put the one before the 

# derive Birthday from IC file
birth_temp = fread("../PDcohorts/PDBP/query_result_InformedConsent_2016-10-26T15-05-40.csv")
birth =   birth_temp[, c(6, 9, 11)]
names(birth) = c("ID", "Age_M", "DATE")
birth = birth %>% 
  mutate(DATE = as.numeric(as.Date(format(DATE, format = "%Y-%m-%d")))) %>% 
  mutate(BRTHDT = DATE - as.numeric(Age_M)*30) %>% 
  select(ID, BRTHDT) %>% distinct(ID, .keep_all = T)

combine2 = left_join(combine1, birth, by="ID") %>% 
  mutate(Age_M = NULL, 
         AGE = (VISIT - BRTHDT)/365.25)

### demographics
## no age at onset
demo_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_Demographics_2018-02-12T15-37-41.csv", header = T)
names(demo_temp)
demo <- demo_temp[, c(6, 5, 4, 10, 12, 13)]

colnames(demo) <- c("ID", "VISIT", "EVENT","FEMALE", "EUROPEAN", "YEARSEDUC1")
demo = demo %>% 
  mutate(FEMALE =  ifelse(FEMALE == "Female", 1, 0),
         EUROPEAN = ifelse(EUROPEAN == "Caucasian (e.g., British Isles, Germany, Peninsular Spain, Latin America, France, Italy, Ireland, Sweden, etc.", 1,0))
demo$YEARSEDUC =  8
demo$YEARSEDUC[demo$YEARSEDUC1 == "9th Grade"] <- 9
demo$YEARSEDUC[demo$YEARSEDUC1 == "11th Grade"] <- 11
demo$YEARSEDUC[demo$YEARSEDUC1 == "Associate degree: academic program"] <- 14
demo$YEARSEDUC[demo$YEARSEDUC1 == "Bachelor's degree (e.g., BA, AB, BS, BBA)"] <- 16
demo$YEARSEDUC[demo$YEARSEDUC1 == "Doctoral degree (e.g., PhD, EdD)"] <- 20
demo$YEARSEDUC[demo$YEARSEDUC1 == "GED or equivalent"] <- 10
demo$YEARSEDUC[demo$YEARSEDUC1 == "High school graduate"] <- 12
demo$YEARSEDUC[demo$YEARSEDUC1 == "Master's degree (e.g., MA, MS, MEng, MEd, MBA)"] <- 22
demo$YEARSEDUC[demo$YEARSEDUC1 == "Some college, no degree"] <- 14

demoI = usuful_func(demo)
# Take the most filled data
demoII = demoI %>% 
  arrange(MON) %>% 
  group_by(ID) %>% 
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>% 
  data.frame %>% 
  arrange(desc(MON)) %>% 
  distinct(ID, .keep_all = T) %>% select(ID, FEMALE, EUROPEAN, YEARSEDUC)

combine3 = left_join(combine2, demoII, by = "ID") 

## age at onset
neuro_temp <- read.csv("../PDcohorts/PDBP/PDBP_data_2018/query_result_NeurologicalExam_2018-03-28T11-45-57.csv", header = T)
neuro <- neuro_temp[, c(6, 5, 4, 11, 12, 13)]
colnames(neuro) <- c("ID", "VISIT", "EVENT", "Dx", "AAO.y", "AAO.m")
neuro = neuro %>% filter(AAO.y>20) %>% 
  mutate(AAO = ifelse(is.na(AAO.m), AAO.y, AAO.y + AAO.m/12)) %>% 
  usuful_func(.) %>% 
  arrange(MON) %>% 
  group_by(ID) %>% 
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>% 
  data.frame %>% 
  arrange(desc(MON)) %>% 
  distinct(ID, .keep_all = T) %>% select(ID, Dx, AAO)

combine4 = left_join(combine3, neuro, by="ID") 

# Family_Hisoty 1st degree
fh_temp = fread("../PDcohorts/PDBP/query_result_FamilyHistory_2016-10-26T15-05-18.csv", na.strings = "")
fh = fh_temp[, c(6,5,4,64:66)]
names(fh)=c("ID", "VISIT", "EVENT", "TYPE", "IND", "RELTVTYP")
# words list
fhtemp = fh %>% filter(!is.na(RELTVTYP)) %>% 
  distinct(RELTVTYP) %>% t %>% as.vector %>% 
  strsplit(., split=",") %>% unlist %>%
  strsplit(., split=";") %>% unlist %>% 
  .[!duplicated(.)]
FH1 = fhtemp %>% # first degree relative
  .[!grepl("great|grand|self|uncle|aunt|wife|cousin|niece", ., ignore.case = T)]
for (i in 1:nrow(fh)){
  x = fh$RELTVTYP[i]
  x_FH = unlist(strsplit(unlist(strsplit(x, split=",")), split = ";"))
  if(sum(x_FH %in% FH1, na.rm = T)>0){fh$FAMILY_HISTORY[i] =1}
  else(fh$FAMILY_HISTORY[i]=0)
}
fh = fh %>% filter(!is.na(ID)) %>% 
  usuful_func(.) %>% 
  select(ID, FAMILY_HISTORY)

combine5 = left_join(combine4, fh, by="ID") %>% 
  mutate(FAMILY_HISTORY = ifelse(is.na(FAMILY_HISTORY), 0, FAMILY_HISTORY))

final = combine5 %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(EVENT=="Baseline") %>% 
  select(ID) %>% # the only ones that the first visit date is the baseline
  inner_join(combine5, ., by ="ID") %>% 
  arrange(AGE) %>% 
  group_by(ID) %>% 
  mutate(AGEatBL = first(AGE)) %>% 
  data.frame

lt = final %>% filter(Dx =="Parkinson's disease") %>% 
  mutate(STUDY_NAME = "PDBP",
         HY3 = ifelse(HY>=3, 1, 0),
         HT=NA, 
         MDS_UPDRS = rowSums(.[, grep("UPDRS", names(.))], na.rm=T),
         oldUPDRS=NA, 
         MMSE=NA, 
         DATE=VISIT, 
         BLDfDIAG = AGEatBL - AAO, 
         STUDY_NAME = "PDBP",
         AD = NA) %>% 
  arrange(DATE) %>% 
  group_by(ID) %>% 
  mutate(TSTART = DATE - first(DATE)) %>% 
  data.frame %>% #DAys from the first visit
  distinct(ID, TSTART, .keep_all=T)

Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (MDS_UPDRS - mean(Baseline$MDS_UPDRS, na.rm=T))/sd(Baseline$MDS_UPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4))) %>% 
  select(VARS)

PDBP = lt %>% select(VARS)


# PICNICS
############################################################
lt = fread("../PDcohorts/picinics/picinics_longitudinal.csv", header = T) %>% 
  mutate(DATE=as.Date(VISIT, format="%m/%d/%Y")) %>% 
  rename(UPDRS1 = UPDRS_Subscale_1, 
         UPDRS2 = UPDRS_Subscale_2, 
         UPDRS3 = UPDRS_Subscale_3,
         UPDRS4 = UPDRS_Subscale_4) %>% 
  mutate(RBD=NA,
         oldUPDRS = NA, 
         HY3 = ifelse(HY>=3.0, 1, 0)) %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame


# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (0 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (MDS_UPDRS - mean(Baseline$MDS_UPDRS, na.rm=T))/sd(Baseline$MDS_UPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4)))

cs =  fread("../PDcohorts/picinics/picinics_crosssectional.csv", header = T) %>%
  .[,1:11] %>% 
  mutate(AGEatBL = AGE,
         BLDfDIAG = AGE-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY=='Caucasian', 1, 0)) %>% 
  distinct(ID, .keep_all=T)

PICNICS = inner_join(lt, cs, by = 'ID') %>% mutate(STUDY_NAME = "PICNICS") %>% select(VARS)

# PPMI
#########################################################
# Data-processing steps should be done prior.
# datastep01, datastep 02, datastep04 is required.
data_folder = "../PDcohorts/PPMI/PPMI180314"
# Basic Characteristics
identity = fread("../PDcohorts/PPMI/outputs/02_LatestDiagDemog.csv",colClasses = c("ID"="character")) %>% 
  rename(PATNO = ID) %>% 
  data.frame %>% 
  .[grep("iPD_YsImg_PPMI", .$DIAG), ] %>% 
  filter(!is.na(DIAGDATE)) %>% 
  select(PATNO, FEMALE, BIRTHDT, EUROPEAN, DIAGDATE, YEARSEDUC, FAMILY_HISTORY)

IDs = identity %>% select(PATNO)

# EVENT_DATE reference
ref = fread("../PDcohorts/PPMI/outputs/01_extract_EVENTDATE.csv",colClasses = c("PATNO"="character")) %>% 
  select(PATNO, EVENT_ID, DATE, Vord) %>% 
  mutate(DATE = as.numeric(as.Date(DATE))) %>% 
  rename(V = Vord)

# LONGITUDINAL DATA
# Combine BL and SC to make a BL
setBL_func = function(data){
  BL = data %>% filter(V <= 0) %>%
    arrange(V) %>% 
    group_by(PATNO) %>% 
    #LOCF for NAs in the UPDRS items for BL if SC is available
    mutate_all(funs(na.locf(., na.rm=F))) %>% 
    data.frame %>% 
    arrange(desc(V)) %>%  
    distinct(PATNO, .keep_all = T) %>% # Keep the latest date among BL or SC
    data.frame %>% 
    mutate(V=-0.1) # Baseline = 0.01
  
  data %>% filter(V>0) %>% 
    rbind(BL, .) %>% 
    mutate(Vord = V, V=NULL)
}


# UPDRS for PD, HY,  and medication data
# Import the output form note2
UPDRS = fread("../PDcohorts/PPMI/outputs/04_UPDRS.csv",
              colClasses = c("PATNO"="character")) %>% 
  # UPDRS4 is often missing at baseline because they don't have motor fluctuations. _> replace 0
  mutate_at(vars(grep("NP4", names(.))), funs(ifelse(is.na(.) & (Vord == -1 | Vord == 0), 0, .))) %>% 
  mutate(UPDRS1 = rowSums(.[,grep('NP1', names(.))], na.rm=F),
         UPDRS2 = rowSums(.[,grep('NP2', names(.))], na.rm=F),
         UPDRS3 = rowSums(.[,grep('NP3', names(.))], na.rm=F),
         UPDRS4 = rowSums(.[,grep('NP4', names(.))], na.rm=F),
         CONST = if_else(NP1CNST>0, 1, 0),
         INS = if_else(NP1SLPN>1, 1, 0),
         MOTORFLUX = if_else(NP4OFF>0, 1, 0), 
         DYSKINESIAS = if_else(NP4WDYSK>0, 1, 0),
         HY3 = ifelse(HY>=3, 1, 0)) %>%
  select(PATNO, DATE, V, UPDRS1, UPDRS2, UPDRS3, UPDRS4, HY, HY3, CONST, INS, MOTORFLUX, DYSKINESIAS) %>% 
  arrange(DATE) %>% 
  group_by(PATNO) %>% 
  mutate_at(vars(names(.)[-1:-6]), funs(na.locf(., na.rm = F))) %>% data.frame %>% 
  mutate(MDS_UPDRS = rowSums(.[,grep('UPDRS', names(.))], na.rm=F))


# New read_func which also returns the missing Date
read_func <- function(data_folder, file){
  data = fread(paste(data_folder, file, sep='/'), na.strings = '', header = T, quote="\"", 
               colClasses = c("PATNO"="character")) %>% 
    inner_join(., IDs, by="PATNO") %>% 
    left_join(., ref, by=c("PATNO", "EVENT_ID")) %>% # give DATE
    arrange(PATNO, DATE)
}

moca=read_func(data_folder, 'Montreal_Cognitive_Assessment__MoCA_.csv') %>% 
  mutate(MOCA = as.numeric(MCATOT)) %>%
  mutate(DEMENTIA = if_else(MOCA < 24, 1, 0)) %>% 
  select(PATNO, V, DATE, MOCA, DEMENTIA)

#SEADL
seadl=read_func(data_folder, 'Modified_Schwab_+_England_ADL.csv') %>% 
  mutate(SEADL = as.numeric(MSEADLG)) %>% 
  select(PATNO, V, DATE, SEADL) 

# HYPOSMIA_longitudinal
upsit=read_func(data_folder, 'University_of_Pennsylvania_Smell_ID_Test.csv') %>% 
  mutate_at(vars(starts_with('UPSITBK')), funs(as.numeric)) %>% 
  mutate(UPSIT = rowSums(.[grep("UPSITBK", names(.))], na.rm = T)) %>% 
  mutate(HYPOSMIA = if_else(UPSIT<21, 1, 0)) %>% 
  select(PATNO, V, DATE, HYPOSMIA)

# DEPR: short GDS > 4 PMC2040268
depr=read_func(data_folder, 'Geriatric_Depression_Scale__Short_.csv') %>% 
  mutate_at(vars(starts_with('GDS')), funs(as.numeric)) %>% 
  mutate(GDS_total =rowSums(.[,grep('GDS*', colnames(.))], na.rm = T)) %>% 
  mutate(DEPR = if_else(GDS_total>5, 1, 0)) %>%
  select(PATNO, V, DATE, DEPR) 

# RLS RBDSQ: RLS(==1), RBD (RBDSQ>5)
rls=read_func(data_folder, 'REM_Sleep_Disorder_Questionnaire.csv') %>%
  mutate_at(c(8:27), funs(as.numeric)) %>% # RBSQ others won't be used.
  mutate(RL = if_else(RLS=='1', 1, 0), 
         q10= if_else(rowSums(.[20:27], na.rm = T)>0, 1, 0)) %>% #comorbidity
  mutate(RBDSQ = rowSums(.[8:19], na.rm = T) + q10) %>% 
  mutate(RBD = if_else(RBDSQ>5, 1, 0)) %>% 
  select(PATNO, V, DATE, RL, RBD)

# Daytime sleepiness (EPworth, score>=10)
sleep=read_func(data_folder, 'Epworth_Sleepiness_Scale.csv') %>%
  mutate_at(vars(starts_with('ESS')), funs(as.numeric)) %>% 
  mutate(SLEEP = if_else(rowSums(.[grep('ESS.$', colnames(.))], na.rm = T)> 9, 1, 0)) %>% 
  select(PATNO, V, DATE, SLEEP) 

# Levodopa or Dopamin Agonist
drug=read_func(data_folder, 'Use_of_PD_Medication.csv') %>%
  mutate(DOPA = if_else(ONLDOPA=='1', 1, 0),
         AGONIST = if_else(ONDOPAG == '1', 1, 0)) %>% 
  select(PATNO, V, DATE, DOPA, AGONIST)


# Longitudinal dataset ##########################################
ltALL = full_join(UPDRS, depr, by=c('PATNO', 'V', "DATE")) %>%
  full_join(., moca, by=c('PATNO', 'V', "DATE")) %>%
  full_join(., rls, by=c('PATNO', 'V', "DATE")) %>%
  full_join(., seadl, by=c('PATNO', 'V', "DATE")) %>%
  full_join(., sleep, by=c('PATNO', 'V', "DATE")) %>%
  full_join(., upsit, by=c('PATNO', 'V', "DATE")) %>% 
  full_join(., drug, by=c("PATNO", "V", "DATE")) %>% 
  # DOPA and AGONIST should be 0 when NA. 
  mutate_at(vars(c('DOPA', 'AGONIST')),funs(if_else(is.na(.), 0, .))) %>%   
  setBL_func(.) %>% 
  arrange(Vord) %>% 
  group_by(PATNO) %>% 
  mutate(DOPA = na.locf(DOPA),
         AGONIST = na.locf(AGONIST)) %>% 
  data.frame

lt = inner_join(identity, ltALL, by="PATNO") %>% 
  mutate(AAO = (DIAGDATE - BIRTHDT)/365.25, 
         AGE = (DATE - BIRTHDT)/365.25, 
         UPDRS4_scaled = as.numeric(scale(UPDRS4)),
         # Unavailable columns
         AD = NA, HT=NA, oldUPDRS=NA, UPDRS_scaled=NA, UPDRS3_scaled=NA, MMSE=NA, PDQ39=NA) %>%
  # Interval from the last observation
  arrange(Vord) %>% 
  group_by(PATNO) %>% 
  mutate(ID = PATNO,
         AGEatBL = first(AGE),
         BLDfDIAG = as.numeric(first(DATE) - DIAGDATE)/365.25, # Baseline Date from DIAGNOSIS (YEARS)
         TSTART = as.numeric(DATE-first(DATE))) %>%  #Visit date from baseline date (DAYS) 
  data.frame

Baseline = subset(lt, Vord == -0.1)

lt = lt %>% 
  mutate(UPDRS_scaled = (MDS_UPDRS - mean(Baseline$MDS_UPDRS, na.rm=T))/sd(Baseline$MDS_UPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         STUDY_NAME = 'PPMI') 

write.csv(lt, paste(data_folder,"_",Sys.Date(), ".csv", sep = ""), row.names = F)

lt =  fread(paste(data_folder,"_",Sys.Date(), ".csv", sep = ""), header = T) %>% data.frame %>% 
  mutate(ID = as.numeric(PATNO)) %>% 
  group_by(ID, TSTART) %>% mutate(dup = n()) %>% 
  data.frame
lt %>% filter(dup>1) %>%
  select(ID,Vord, TSTART) %>% data.frame %>% arrange(ID)
# For the observation in the same day (month), use na.locf

PPMI = lt %>% filter(dup>1) %>% 
  group_by(ID, TSTART) %>% mutate_all(.,funs(na.locf(., na.rm = F))) %>% 
  data.frame %>% 
  distinct(ID, TSTART, .keep_all = T) %>% 
  rbind(., subset(lt, dup==1)) %>% 
  mutate(dup=NULL) %>% 
  arrange(ID, TSTART) %>% 
  select(VARS)




# PRECEPT 
########################################################
lt = fread("../PDcohorts/PRECEPT/precept_longitudinal.csv", header = T) %>% 
  mutate(DATE=as.Date(VISIT, format="%d-%b-%y")) %>% 
  rename(UPDRS1 = updrsmental, 
         UPDRS2 = updrsadl,
         UPDRS3 = updrsmotor) %>% 
  mutate(oldUPDRS = rowSums(.[,grep('UPDRS', colnames(.))], na.rm = T),
         HYPOSMIA = if_else(UPSIT < 21, 1, 0),
         RL = NA,
         HT = NA, 
         RBD = NA,
         UPDRS4 = NA,
         MDS_UPDRS = NA,
         UPDRS4=NA, 
         HY3 = ifelse(HY>=3.0, 1, 0),
         DEMENTIA = ifelse(MMSE<27, 1, 0)) %>% # MOCA is also available but more missing. 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame

# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (0 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (oldUPDRS - mean(Baseline$oldUPDRS, na.rm=T))/sd(Baseline$oldUPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4))) 

cs =  fread("../PDcohorts/PRECEPT/precept_crosssectional.csv", header = T) %>%
  mutate(AGEatBL = Age,
         BLDfDIAG = Age-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY=='Caucasian', 1, 0)) %>% 
  mutate(STUDY_NAME = 'PRECEPT') 

PRECEPT = inner_join(lt, cs, by = 'ID') %>% select(VARS)

# Predict
######################################################################
# The cohort was removed.
# Predict = fread('../PDcohorts/Predict/data_on_PREDICT_subjects_sent_to_PPP.csv') %>% na_if(999)
# 
# Predict_basic = Predict %>% 
#   mutate(SUBJID = id,
#          FEMALE = 1-gender_bin,
#          AGEatBL = age_0,
#          FAMILY_HISTORY = efirstdegree_0) %>% 
#   select(SUBJID, FEMALE, AGEatBL, FAMILY_HISTORY)
# 
# Predict_f0 = Predict %>% 
#   mutate(SUBJID=id, 
#          Case = 0,
#          Age = age_0,
#          HT = eheadinjury_0,
#          DEPR = ehadsmoderate_0,
#          CONST = econst_0,
#          UPSIT = upsit_0,
#          RBD = rbd_0,         
#          MOCA = moca, 
#          DEMENTIA = ifelse(moca<24, 1, 0),
#          UPDRS_scaled = as.numeric(scale(blindupdrs))) %>% # I don't know which UPDRS it is.
#   select(SUBJID, Case, Age, HT, DEPR, CONST, UPSIT, RBD, MOCA, DEMENTIA, UPDRS_scaled)
# 
# Predict_f1 = Predict %>% 
#   mutate(SUBJID=id, 
#          Case = pd_1,
#          Age = age_1,
#          HT = eheadinjury_1,
#          DEPR = ehadsmoderate_1,
#          CONST = econst_1,
#          UPSIT = NA,
#          RBD = rbd_1, 
#          MOCA=NA, DEMENTIA=NA, UPDRS_scaled=NA) %>% 
#   select(SUBJID, Case, Age, HT, DEPR, CONST, UPSIT, RBD, MOCA, DEMENTIA, UPDRS_scaled)
# 
# Predict_f2 = Predict %>% 
#   mutate(SUBJID=id, 
#          Case = pd_2,
#          Age = age_2,
#          HT = eheadinjury_2,
#          DEPR = ehadsmoderate_2,
#          CONST = econst_2,
#          UPSIT = NA,
#          RBD = rbd_2, 
#          MOCA=NA, DEMENTIA=NA, UPDRS_scaled=NA) %>% 
#   select(SUBJID, Case, Age, HT, DEPR, CONST, UPSIT, RBD, MOCA, DEMENTIA, UPDRS_scaled)
# 
# Predict_f3 = Predict %>% 
#   mutate(SUBJID=id, 
#          Case = pd_3,
#          Age = age_3,
#          HT = eheadinjury_3,
#          DEPR = ehadsmoderate_3,
#          CONST = econst_3,
#          UPSIT = upsit_3,
#          RBD = rbd_3, 
#          MOCA=NA, DEMENTIA=NA, UPDRS_scaled=NA) %>% 
#   select(SUBJID, Case, Age, HT, DEPR, CONST, UPSIT, RBD, MOCA, DEMENTIA, UPDRS_scaled)
# 
# PREDICT = bind_rows(Predict_f0, Predict_f1, Predict_f2, Predict_f3) %>% 
#   left_join(., Predict_basic, by = "SUBJID") %>% 
#   group_by(SUBJID) %>% 
#   arrange(Age) %>% 
#   mutate(TSTART = floor((Age - first(Age))*365.25)) %>% 
#   data.frame %>% 
#   mutate(STUDY_NAME = "PREDICT",
#          HYPOSMIA = ifelse(UPSIT<=21, 1, 0), 
#          RBD = ifelse(RBD>=5, 1, 0), 
#          MOTORFLUX = NA, DYSKINESIAS = NA, RL = NA, SLEEP = NA, 
#          INS = NA, HY3 = NA, EUROPEAN = NA, DOPA = NA, AGONIST = NA, 
#          AAO=NA, AD=NA, YEARSEDUC=NA, BLDfDIAG=NA,
#          HY=NA,UPDRS1_scaled=NA,UPDRS2_scaled=NA,UPDRS3_scaled=NA,UPDRS4_scaled=NA,
#          UPDRS1=NA, UPDRS2=NA, UPDRS3=NA, UPDRS4=NA, MDS_UPDRS=NA, oldUPDRS=NA, 
#          MMSE=NA,SEADL=NA,PDQ39=NA) %>% 
#   rename(ID = SUBJID) %>% 
#   select(VARS)



# SCOPA
#######################################################################
lt = fread("../PDcohorts/SCOPA/SCOPA_longitudinal.csv", header = T) %>% 
  mutate(DATE=as.Date(VISIT, format="%m/%d/%Y")) %>% 
  rename(UPDRS1 = UPDRS_Subscale_1, 
         UPDRS2 = UPDRS_Subscale_2, 
         UPDRS3 = UPDRS_Subscale_3,
         UPDRS4 = UPDRS_Subscale_4) %>% 
  mutate(RBD=NA,
         oldUPDRS = NA, 
         HY3 = ifelse(HY>=3.0, 1, 0)) %>%
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame


# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (0 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

lt = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() == 1) %>% 
  mutate(idx=1) %>% data.frame %>% #idx == number of the observation in the same ID, Date
  rbind(., lt_process_duplicate) %>% 
  mutate(UPDRS4 = as.numeric(scale(UPDRS4)))




###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (MDS_UPDRS - mean(Baseline$MDS_UPDRS, na.rm=T))/sd(Baseline$MDS_UPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4))) 


cs =  fread("../PDcohorts/SCOPA/SCOPA_crosssectional.csv", header = T) %>%
  .[,1:11] %>% 
  mutate(AGEatBL = AGE,
         BLDfDIAG = AGE-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY=='Caucasian', 1, 0)) %>% 
  distinct(ID, .keep_all=T)

SCOPA = inner_join(lt, cs, by = 'ID') %>% mutate(STUDY_NAME = "SCOPA") %>% select(VARS)


# UDALL
############################################################################
lt = fread("../PDcohorts/UDALL/penn_longitudinal.csv", header = T) %>% 
  mutate(DATE = as.Date(VISIT, '%m/%d/%Y')) %>% 
  rename(UPDRS1 = UPDRS_Subscale_1, 
         UPDRS2 = UPDRS_Subscale_2, 
         UPDRS3 = UPDRS_Subscale_3,
         UPDRS4 = UPDRS_Subscale_4,
         MOCA = MoCA,
         DEMENTIA = DEMENT,
         DEPR = DEPRESS) %>% 
  mutate(oldUPDRS = rowSums(.[grep('UPDRS',colnames(.))], na.rm = T),
         RBD=NA,
         HYPOSMIA=NA, # will be replaced later
         DOPA = ifelse(DOPA>0, 1, DOPA),
         AGONIST = ifelse(AGONIST>0, 1, AGONIST), 
         HY3 = ifelse(HY>=3.0, 1, 0)) %>% 
  group_by(ID) %>% 
  arrange(DATE) %>% 
  mutate(TSTART = as.numeric(DATE - first(DATE))) %>%  #DAys from the first visit
  data.frame

# Duplicated observation were kept using LOCF-like method.
lt_process_duplicate = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() > 1) %>%  # Filter only duplicated ones (0 case)
  mutate_all(funs(na.locf(.,na.rm = FALSE))) %>%
  mutate(idx = row_number()) %>% #Index within group
  filter(idx==max(idx)) %>% #Keep the last one
  data.frame

lt = lt %>% 
  group_by(ID, DATE) %>% 
  filter(n() == 1) %>% 
  mutate(idx=1) %>% data.frame %>% #idx == number of the observation in the same ID, Date
  rbind(., lt_process_duplicate)


###Construct Baseline data and longitudinal data####
Baseline = lt %>% 
  filter(TSTART == 0)

lt = lt %>% 
  mutate(UPDRS_scaled = (oldUPDRS - mean(Baseline$oldUPDRS, na.rm=T))/sd(Baseline$oldUPDRS, na.rm = T),
         UPDRS1_scaled = (UPDRS1 - mean(Baseline$UPDRS1, na.rm=T))/sd(Baseline$UPDRS1, na.rm = T),
         UPDRS2_scaled = (UPDRS2 - mean(Baseline$UPDRS2, na.rm=T))/sd(Baseline$UPDRS2, na.rm = T),
         UPDRS3_scaled = (UPDRS3 - mean(Baseline$UPDRS3, na.rm=T))/sd(Baseline$UPDRS3, na.rm = T),
         UPDRS4_scaled = as.numeric(scale(UPDRS4)),
         HT = ifelse(HT==9, NA, HT)) 


cs =  fread("../PDcohorts/UDALL/penn_crosssectional.csv", header = T) %>% .[,1:11] %>% 
  mutate(AGEatBL = AGE,
         BLDfDIAG = AGE-AAO, # Years from diagnosis to baseline (Not available in this dataset)
         EUROPEAN = if_else(ANCESTRY== 'European', 1, 0)) %>% 
  distinct(ID, .keep_all=T) %>% 
  left_join(.,fread("../PDcohorts/UDALL/Cross sectional UPSIT data.csv"), by = c("ID"="INDDID")) %>% 
  mutate(HYPOSMIA = ifelse(`UPSIT score` <21, 1, 0)) %>% 
  mutate(`UPSIT score`=NULL) %>% 
  distinct(ID, .keep_all = T) # UPSIT has two obs for INDDID 114740

UDALL = cs %>% select(-HYPOSMIA) %>% inner_join(., lt, by = 'ID') %>% mutate(STUDY_NAME = "PRECEPT") %>%  select(VARS)





###############################################################
#################################################################




Cohort13 = rbind(CORIELL, DATATOP, DIGPD, HBS, OSLO_LT, PARKFIT, PARKWEST, PDBP, PICNICS, PPMI, PRECEPT, SCOPA, UDALL) %>% 
  arrange(STUDY_NAME, ID, TSTART) %>% 
  mutate(YEARfDIAG = BLDfDIAG+TSTART/365.25)

write.csv(Cohort13, paste("Pheno13_", Sys.Date(), ".csv"), row.names = F)

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
