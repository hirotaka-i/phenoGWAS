# Create phenofile for this analysis genetics data
# Original phenofile created from raw cohort data and only genotyped individuals were selected.
# HBS data has not genetic-phenotypic table so cannot join here.
rm(list=ls())
require(openxlsx)
require(dplyr)
require(data.table)
require(zoo)
###################################
# Step1
# Read IDs for genotyped individuals
# Add phenotype ID column to match the phenotype file
# Sheet 1 all except Oslo individuals
# Sheet 2 Oslo individuals
###################################
temp1 = read.xlsx("data/progression_GWAS_imputed_IDs.xlsx", sheet = 1) %>% 
  arrange(DATASET, ID)
temp2 = temp1 %>% 
  group_by(DATASET, ID) %>%
  mutate(PID = case_when(
    DATASET != "HBS"~ strsplit(ID, "_")[[1]][1],
    strsplit(ID, "_")[[1]][2] == "PD" ~ strsplit(ID, "_")[[1]][3],
    TRUE~ strsplit(ID, "_")[[1]][2])) %>% data.frame %>%
  # Need to convert the ID for combining with fam/bin/bed files
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

print("complete step1")
#############################
# Step2
# Match above file with phenotype file
# PhenoFile should be in PDcohorts/codes. Before retrieve, make sure the file is the latest
#############################
temp1 = fread("outputs/ID_imputed.csv")
Pheno = fread("data/Pheno13_2018-06-29.csv", stringsAsFactors = F) %>% 
  rename(COHORT = STUDY_NAME) %>% 
  mutate(ID = case_when(
    COHORT == "PARKFIT"~ tolower(ID), 
    TRUE~ID
  ))

PhenoFile = inner_join(temp1, Pheno, by=c("COHORT", "PID"="ID"))
write.csv(PhenoFile, "outputs/PhenoFile.csv", row.names = F)

rm(list=ls())
print("compete step2")

#############################
# Step3 
## prepare for rvtest (baseline lenear regression (lm) / logistic regression(glm))
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
ltg_all=fread('outputs/PhenoFile.csv') %>% arrange(COHORT, ID, TSTART) %>% data.frame
COHORTs = unique(ltg_all$COHORT)
COHORTs
names(ltg_all)
BinomT = names(ltg_all)[7:17]
BinomT
ContT = names(ltg_all)[18:32]
ContT
ContT = ContT[-(7:12)] # erase UPDRS1-4, oldUPDRS, MDSUPDRS
ContT
COVs = c('FEMALE', 'YEARSEDUC', "FAMILY_HISTORY", 'AAO', "BLDfDIAG",  'DOPA', 'AGONIST')

plink_ltg = ltg_all %>% 
  mutate(FID = ID,
         IID = ID,
         FATID = 0,
         MATID = 0,
         SEX = 1+FEMALE) # In plink format, SEX 0 missing, 1 male, 2 female

MODELs = data.frame(COHORT_covs_lm_glm_surv = rep(NA, length(COHORTs)))
unlink("outputs/rvtest", recursive = T)
dir.create("outputs/rvtest", showWarnings = F)
unlink("outputs/surv", recursive = T)
dir.create("outputs/surv", showWarnings = F)
unlink("outputs/long", recursive = T)
dir.create("outputs/long", showWarnings = F)

PCAs = "/data/LNG/Hirotaka/progGWAS/pca/"

for(i in 1:length(COHORTs)){
  BL = plink_ltg %>% filter(TSTART ==0 & COHORT == COHORTs[i]) %>% 
    mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1))) # 0 missing, 1 control, 2 case
  if(COHORTs[i]=="OSLO"){ # Oslo has no appearent baseline)
    BL = plink_ltg %>% filter (COHORT == COHORTs[i]) %>% filter(-365 < TSTART & TSTART < 365) %>% 
      distinct(ID, .keep_all = T) %>% 
      mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1)))
  }
  # PC1-5
  PC = fread(paste(PCAs, COHORTs[i], ".pca.eigenvec", sep = "")) %>% select(-FID) %>% distinct(IID, .keep_all = T) %>% mutate(IID = as.character(IID))
  # Binomial outcomes
  binom = BL %>% select(FID, IID, FATID, MATID, SEX, BinomT)
  # Avaial phenos for cohort_i 
  temp = binom[,6:ncol(binom)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  BinomTs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  if(length(BinomTs_i) != 0){ # if there are such variables..
    binom_new = binom %>% select(FID, IID, FATID, MATID, SEX, BinomTs_i)
    write.table(binom_new, paste("outputs/rvtest/", COHORTs[i], ".binom", sep = ""), row.names = F, quote = F, sep = "\t")
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
  write.table(cov_new, paste("outputs/rvtest/", COHORTs[i], ".cov", sep = ""), row.names = F, quote = F, sep = "\t")
  
  # Continuous outcomes (use rvtest only when longitudinal data is not available)
  cont = plink_ltg %>% filter(COHORT == COHORTs[i])
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
    write.table(cont_new, paste("outputs/rvtest/", COHORTs[i], ".contbl", sep = ""), row.names = F, quote = F, sep = "\t")
  }else{ContTs_i_cs = ""}
    # Longitudinal set: it also should have covsariates because they will be analyzed with outcomes in R
  if(length(ContTs_i_lt)!=0){
    COVs_i_lt_temp = gsub("BLDfDIAG", "YEARfDIAG", COVs)# replace BLDfDIAG -> YEARfDIAG
    # Filger out non-effective COVs_i_lt
    temp = cont %>% select(COVs_i_lt_temp) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
    temp$NAMES = rownames(temp) # rownames were converted to the covariates.
    COVs_i_lt = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector
    cont_new = cont %>% select(FID, IID, FATID, MATID, SEX, COVs_i_lt, ContTs_i_lt) %>% left_join(., PC, by = "IID")
    write.table(cont_new, paste("outputs/long/", COHORTs[i], ".cont", sep = ""), row.names = F, quote = F, sep = "\t")
  }else{ContTs_i_lt = ""}
  # Survival outcome
  PHENOSURV = c()
  for (j in (1:length(BinomT))){
    ltg = ltg_all %>% filter(COHORT == COHORTs[i])
    ltg$OUTCOME = ltg[,BinomT[j]]
    
    cov2in = ltg %>% select(COVs) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
    cov2in$NAMES = rownames(cov2in) # rownames were converted to the covariates.
    COVs2in =  cov2in %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
    
    ltg = ltg %>% select(ID, TSTART, OUTCOME, COVs2in)
    time2follow = ltg %>%
      arrange(TSTART) %>%
      filter(!is.na(OUTCOME)) %>%
      group_by(ID, .keep_all = T) %>%
      mutate(BEGIN = first(TSTART),
             END = last(TSTART)) %>% # BEGIN and END are identical within each ID
      mutate(END = ifelse(OUTCOME==1, TSTART, END)) %>% # If 1, then END will be updated
      data.frame %>%
      arrange(END) %>% 
      distinct(ID, .keep_all = T) %>%  # take the smallest END
      select(ID, BEGIN, END)
    cohort_by_first_event = left_join(ltg, time2follow, by = 'ID') %>%
      filter(BEGIN <= TSTART & TSTART <= END) %>% # filter observation between BEGIN and END
      mutate(OUTCOME = ifelse(is.na(OUTCOME), 0, OUTCOME)) # NAs in the period replaced by 0
    
    cohort = cohort_by_first_event %>% 
      arrange(TSTART) %>% 
      group_by(ID) %>%
      mutate(BLDfDIAG = BLDfDIAG + TSTART/365.25, # Follow-up starts TSTART > 0 sometimes, Needs update
             TSTART = TSTART - first(TSTART)) %>% # New baseline set 0 (Especially, ParkWest)
      mutate(TSTOP = lead(TSTART), # The TSTART[i+1] will be the stop time
             OUTCOME = lead(OUTCOME)) %>%  # the OUTCOME is updated by the one at the stop time
      mutate(TSTOP = ifelse(is.na(TSTOP), END, TSTOP)) %>% # replace the 
      data.frame %>% arrange(ID, TSTART) %>% 
      filter(!is.na(OUTCOME)) # the basially last follow-up data will be erased (Right censored)
    
    if(nrow(cohort)!=0 & var(cohort$OUTCOME, na.rm = T) != 0){
      cohort %>% left_join(., PC, by = c("ID" = "IID")) %>% 
        write.table(., paste("outputs/surv/", BinomT[j], ".", COHORTs[i], ".surv", sep = ""), row.names = F, quote = F, sep = "\t")
      PHENOSURV = c(PHENOSURV, BinomT[j])
    }
  }
  MODELs[i,] = paste(COHORTs[i], paste(COVs_i, collapse = ","), paste(ContTs_i_cs, collapse = ","), paste(BinomTs_i, collapse = ","), 
                     paste(PHENOSURV, collapse = ","), paste(ContTs_i_lt, collapse = ","), sep=":")
}
write.table(MODELs, "outputs/MODELs.txt", row.names = F, quote = F)


# ID conversion for analysis with Imputed data
## CORIELL, DATATOP, PARKFIT, PARKWEST, PICNICS, PRECEPT, SCOPA -> ID_ID
## OSLO, PDBP, PPMI: OK
## HBS: HBS_PD_INVDY797MD3  HBS_PD_INVXJ077RRF

for(FILE in list.files(path = "outputs/rvtest")){
  if(grepl("CORIELL|DATATOP|PARKFIT|PARKWEST|PICNICS|PRECEPT|SCOPA", FILE)){
    fread(paste("outputs/rvtest/", FILE, sep = "")) %>%
      mutate(FID = paste(FID, FID, sep="_"),
             IID = paste(IID, IID, sep="_")) %>% 
      write.table(., paste("outputs/rvtest/", FILE, sep = ""), row.names = F, quote = F, sep = "\t")
  }
}

for(FILE in list.files(path = "outputs/surv")){
  if(grepl("CORIELL|DATATOP|PARKFIT|PARKWEST|PICNICS|PRECEPT|SCOPA", FILE)){
    fread(paste("outputs/surv/", FILE, sep = "")) %>%
      mutate(ID = paste(ID, ID, sep="_")) %>% 
      write.table(., paste("outputs/surv/", FILE, sep = ""), row.names = F, quote = F, sep = "\t")
  }
}


print("complete")