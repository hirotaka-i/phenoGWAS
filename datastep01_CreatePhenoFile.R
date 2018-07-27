# Create phenofile for this analysis genetics data
# Original phenofile created from raw cohort data and only genotyped individuals were selected.
# HBS data has not genetic-phenotypic table so cannot join here.
rm(list=ls())
require(dplyr)
require(data.table)
require(zoo)
###################################
# Step1
# Crate IID (genome ID) and ID (phenotype ID) mathcing table
# Default is IID = ID_ID
# OSLO and HBS: genome ID is different from clinical ID
# PARKFIT: genome ID is lower lettered.
# PDBP: genome ID is clinical ID + ST_xxxxx
# DIGPD used two chips. DIGPD_neuroX and DIGPD_chip (two datasets in one cohort)
###################################
ID_g_v1 = fread("_IDlist_imputed.txt", header = F)
names(ID_g_v1) = c("DATASET","IID")
ID_gOK = ID_g_v1 %>% 
  filter(!DATASET %in% c("OSLO", "HBS", "DIGPD_neuroX", "DIGPD_chip")) %>% 
  group_by(DATASET, IID) %>%
  mutate(ID = case_when(
    DATASET=="PARKFIT"~toupper(strsplit(IID, "_")[[1]][2]),
    TRUE~strsplit(IID, "_")[[1]][1])) %>% 
  mutate(COHORT = DATASET) # DIGPD has two datasets in a cohort (later)
ID_gOSLO = fread("data/Oslo_GenotypingClinical_IDs.csv") %>% 
  mutate(DATASET = "OSLO",
         IID = paste(FID, IID, sep = "_")) %>% 
  rename(ID = Clincal_ID) %>% 
  semi_join(., ID_g_v1, by = c("DATASET", "IID")) %>% # keep only imputed
  select(DATASET, IID, ID) %>% 
  mutate(COHORT = DATASET)
ID_gHBS = fread("data/plink.genome") %>% 
  mutate(
    ID = case_when(grepl("HBS", IID1) ~ IID2, TRUE ~ IID1),
    IID= case_when(grepl("HBS", IID1) ~ IID1,TRUE ~ IID2)
  ) %>% select(ID, IID) %>% 
  mutate(DATASET = "HBS") %>% 
  semi_join(., ID_g_v1, by = c("DATASET", "IID")) %>% 
  mutate(COHORT = DATASET)
ID_gDIGPD = ID_g_v1 %>% 
  filter(DATASET %in% c("DIGPD_neuroX", "DIGPD_chip")) %>% 
  group_by(DATASET, IID) %>%
  mutate(ID = strsplit(IID, "_")[[1]][1]) %>% 
  mutate(COHORT = "DIGPD")
  
  

ID_g = bind_rows(ID_gOK, ID_gOSLO, ID_gHBS, ID_gDIGPD)

  
# read phenotype file
ID_p_v1 = fread("data/Pheno13_2018-06-29.csv", stringsAsFactors = F) %>% 
  distinct(STUDY_NAME, ID) %>%
  rename(COHORT = STUDY_NAME)
ID_pg = left_join(ID_g, ID_p_v1, by = c("COHORT", "ID"))
ID_pg %>% filter(is.na(ID)) %>% head # find IID without ID

write.csv(ID_pg, "outputs/ID_imputed.csv", row.names = F)

print("complete step1")




#############################
# Step2
# Match above file with phenotype file
# PhenoFile should be in PDcohorts/codes. Before retrieve, make sure the file is the latest
#############################
temp1 = fread("outputs/ID_imputed.csv")
Pheno = fread("data/Pheno13_2018-06-29.csv", stringsAsFactors = F) %>% 
  rename(COHORT = STUDY_NAME)
PhenoFile = inner_join(temp1, Pheno, by=c("COHORT", "ID")) %>% arrange(DATASET, IID, TSTART)
write.csv(PhenoFile, "outputs/PhenoFile.csv", row.names = F)
rm(list=ls())
print("compete step2")

#############################
# Step3 
## prepare data&model 
### for cross-sectional/baseline analysis) with rvtest
### for logistic analysis  using R

# PCs will be integrated but the file IID is different from IID in imputed file


## prepare for survival test in R
## prepare for linear mixed effect test in abel
### linear regression will be applied only when the data is only cross-sectional.
### Otherwise, lenear mixed model will be applied. 
# .binom: binomial phenotypes file per cohort
# .contbl: continuous phenotype file per cohort (baseline)
# .surv: phenotype and covariates file for survival analysis per cohort
# .cov: covariates file per cohort
# .contlg: continuous phenotype file per cohort (longitudinal)
# .covlg: covariates file per cohort
# MODELs.txt: cohort:covs:outcomes(lm):outcomes(glm):
# MODELs2.txt: cohort:covs:outcomes(survival)
#############################

# Read data
ltg_all=fread('outputs/PhenoFile.csv') %>% arrange(DATASET, ID, TSTART) %>% 
  data.frame %>% 
  mutate(SEADL70 = ifelse(SEADL <= 70, 1, 0))
temp = ltg_all %>% distinct(DATASET, COHORT)
temp
DATASETs = temp$DATASET
COHORTs = temp$COHORT
names(ltg_all)
BinomT = c("CONST","DEMENTIA", "DEPR", "DYSKINESIAS", "HY3", "HYPOSMIA", "INS", "MOTORFLUX", "RBD", "RL", "SLEEP", "SEADL70")
ContT = c("HY", "MMSE", "MOCA", "SEADL", "UPDRS_scaled", "UPDRS1_scaled", "UPDRS2_scaled", "UPDRS3_scaled", "UPDRS4_scaled")
COVs = c('FEMALE', 'YEARSEDUC', "FAMILY_HISTORY", 'AAO', "BLDfDIAG",  'DOPA', 'AGONIST')
plink_ltg = ltg_all %>% 
  mutate(FID = IID,
         IID = IID,
         FATID = 0,
         MATID = 0,
         SEX = 1+FEMALE) # In plink format, SEX 0 missing, 1 male, 2 female

MODELs = data.frame(DATASETs_covs_lm_logistic_surv_lmm = rep(NA, length(DATASETs)))
unlink("outputs/rvtest", recursive = T)
dir.create("outputs/rvtest", showWarnings = F)
unlink("outputs/surv", recursive = T)
dir.create("outputs/surv", showWarnings = F)
unlink("outputs/long", recursive = T)
dir.create("outputs/long", showWarnings = F)

PCAs = "/data/LNG/Hirotaka/progGWAS/pca/"

for(i in 1:length(DATASETs)){
  # Baseline for binomial outcomes
  BL = plink_ltg %>% filter(TSTART ==0 & DATASET == DATASETs[i]) %>% 
    mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1))) # 0 missing, 1 control, 2 case
  if(COHORTs[i]=="OSLO"){ # Oslo has no appearent baseline)
    BL = plink_ltg %>% filter (COHORT == COHORTs[i]) %>% filter(-365 < TSTART & TSTART < 365) %>% 
      distinct(ID, .keep_all = T) %>% 
      mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1)))
  }
  # PC1-5
  PC = fread(paste(PCAs, DATASETs[i], ".pca.eigenvec", sep = "")) 
  names(PC)=c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5")
  PC = PC %>% select(-FID) %>% distinct(IID, .keep_all = T) %>% mutate(IID = as.character(IID))
  # Binomial outcomes
  binom = BL %>% select(FID, IID, FATID, MATID, SEX, BinomT)
  # Avaial phenos for dataset_i 
  temp = binom[,6:ncol(binom)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  BinomTs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  if(length(BinomTs_i) != 0){ # if there are such variables..
    binom_new = binom %>% select(FID, IID, FATID, MATID, SEX, BinomTs_i)
    write.table(binom_new, paste("outputs/rvtest/", DATASETs[i], ".binom", sep = ""), row.names = F, quote = F, sep = "\t")
  }else{
    BinomTs_i_cs = ""
  }
  # Covariate
  cov = BL %>% select(FID, IID, FATID, MATID, SEX, COVs) # Cov can use 0 not as missing.
  # filter out COVARIATES with all the same values or all NAs
  temp = cov[,6:ncol(cov)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  COVs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  cov_new = cov %>% select(FID, IID, FATID, MATID, SEX, COVs_i) %>% left_join(., PC, by = "IID")
  write.table(cov_new, paste("outputs/rvtest/", DATASETs[i], ".cov", sep = ""), row.names = F, quote = F, sep = "\t")
  
  # Continuous outcomes (use rvtest only when longitudinal data is not available)
  cont = plink_ltg %>% filter(DATASET == DATASETs[i])
  # filter out COVARIATES with all the same values or all NAs
  temp = cont %>% select(ContT) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  ContTs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  temp2 = cont %>% filter(TSTART>0) %>% group_by(IID) %>% summarize_at(vars(ContTs_i),funs(sum(!is.na(.)))) %>% 
    summarize_at(vars(ContTs_i),funs(max)) %>% t %>% data.frame # number of observataion per ID (max)
  ContTs_i_cs = setdiff(ContTs_i, rownames(temp2)[temp2[,1] > 0]) # only cross-sectional continuous traits
  ContTs_i_lt=rownames(temp2)[temp2[,1] > 0]
  if(length(ContTs_i_cs) != 0){ # if there are such variables..
    cont_new = cont %>% filter(TSTART == 0) %>% select(FID, IID, FATID, MATID, SEX, ContTs_i_cs) %>% left_join(., PC, by = "IID")
    write.table(cont_new, paste("outputs/rvtest/", DATASETs[i], ".contbl", sep = ""), row.names = F, quote = F, sep = "\t")
  }else{ContTs_i_cs = ""}
  # Longitudinal set: it also should have covsariates because they will be analyzed with outcomes in R
  if(length(ContTs_i_lt)!=0){
    COVs_i_lt_temp = gsub("BLDfDIAG", "YEARfDIAG", COVs)# replace BLDfDIAG -> YEARfDIAG
    # Filger out non-effective COVs_i_lt
    temp = cont %>% select(COVs_i_lt_temp) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
    temp$NAMES = rownames(temp) # rownames were converted to the covariates.
    COVs_i_lt = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector
    cont_new = cont %>% select(FID, IID, FATID, MATID, SEX, COVs_i_lt, ContTs_i_lt) %>% left_join(., PC, by = "IID")
    write.table(cont_new, paste("outputs/long/", DATASETs[i], ".cont", sep = ""), row.names = F, quote = F, sep = "\t")
  }else{ContTs_i_lt = ""}
  # Survival outcome
  PHENOSURV = c()
  for (j in (1:length(BinomT))){
    print(paste("i is ", i , ". And j is", j))
    ltg = ltg_all %>% filter(DATASET == DATASETs[i])
    ltg$OUTCOME = ltg[,BinomT[j]]
    
    cov2in = ltg %>% select(COVs) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
    cov2in$NAMES = rownames(cov2in) # rownames were converted to the covariates.
    COVs2in =  cov2in %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
    
    ltg = ltg %>% select(IID, TSTART, OUTCOME, COVs2in)
    time2follow = ltg %>%
      arrange(TSTART) %>%
      filter(!is.na(OUTCOME)) %>%
      group_by(IID, .keep_all = T) %>%
      mutate(BEGIN = first(TSTART),
             END = last(TSTART)) %>% # BEGIN and END are identical within each ID
      mutate(END = ifelse(OUTCOME==1, TSTART, END)) %>% # If 1, then END will be updated
      data.frame %>%
      arrange(END) %>% 
      distinct(IID, .keep_all = T) %>%  # take the smallest END
      select(IID, BEGIN, END)
    cohort_by_first_event = left_join(ltg, time2follow, by = 'IID') %>%
      filter(BEGIN <= TSTART & TSTART <= END) %>% # filter observation between BEGIN and END
      mutate(OUTCOME = ifelse(is.na(OUTCOME), 0, OUTCOME)) # NAs in the period replaced by 0
    
    cohort = cohort_by_first_event %>% 
      arrange(TSTART) %>% 
      group_by(IID) %>%
      mutate(BLDfDIAG = BLDfDIAG + TSTART/365.25, # Follow-up starts TSTART > 0 sometimes, Needs update
             TSTART = TSTART - first(TSTART)) %>% # New baseline set 0 (Especially, ParkWest)
      mutate(TSTOP = lead(TSTART), # The TSTART[i+1] will be the stop time
             OUTCOME = lead(OUTCOME)) %>%  # the OUTCOME is updated by the one at the stop time
      mutate(TSTOP = ifelse(is.na(TSTOP), END, TSTOP)) %>% # replace the 
      data.frame %>% arrange(IID, TSTART) %>% 
      filter(!is.na(OUTCOME)) # the basially last follow-up data will be erased (Right censored)
    
    if(nrow(cohort)!=0 & var(cohort$OUTCOME, na.rm = T) != 0){
      cohort %>% left_join(., PC, by = c("IID")) %>% 
        write.table(., paste("outputs/surv/", BinomT[j], ".", DATASETs[i], ".surv", sep = ""), row.names = F, quote = F, sep = "\t")
      PHENOSURV = c(PHENOSURV, BinomT[j])
    }
  }
  MODELs[i,] = paste(DATASETs[i], paste(COVs_i, collapse = ","), paste(ContTs_i_cs, collapse = ","), paste(BinomTs_i, collapse = ","), 
                     paste(PHENOSURV, collapse = ","), paste(ContTs_i_lt, collapse = ","), sep=":")
}

write.table(MODELs, "outputs/MODELs.txt", row.names = F, quote = F)
print("all complete")
