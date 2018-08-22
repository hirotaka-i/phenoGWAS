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
ID_g_v1 = fread("t/_IDlist_imputed.txt", header = F)
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
ID_gHBS = fread("data/HBSIDmatch.txt") %>% 
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
cat("find IID without ID. No output is good \n")
ID_pg %>% filter(is.na(ID)) %>% head # find IID without ID

write.csv(ID_pg, "t/ID_imputed.csv", row.names = F)

cat("complete step1 \n")


#############################
# Step2
# Match above file with phenotype file
# PhenoFile should be in PDcohorts/codes. Before retrieve, make sure the file is the latest
#############################
temp1 = fread("t/ID_imputed.csv")
Pheno = fread("data/Pheno13_2018-06-29.csv", stringsAsFactors = F) %>% 
  rename(COHORT = STUDY_NAME)
PhenoFile = inner_join(temp1, Pheno, by=c("COHORT", "ID")) %>% arrange(DATASET, IID, TSTART)
write.csv(PhenoFile, "data/PhenoFile.csv", row.names = F)
rm(list=ls())
cat("compete step2 \n")

#############################
# Step3 
## prepare data&model for analysis 
### for cross/lnsgl use rvtest
### for lncns/lnmxi/lnmxs/coxhm use R
### models.txt is recording a model.
#############################

# Read data
ltg_all=fread('data/PhenoFile.csv') %>% arrange(DATASET, ID, TSTART) %>% 
  data.frame %>% 
  mutate(SEADL70 = ifelse(SEADL <= 70, 1, 0))
temp = ltg_all %>% distinct(DATASET, COHORT)
DATASETs = temp$DATASET
BinomT = c("CONST","DEMENTIA", "DEPR", "DYSKINESIAS", "HY3", "HYPOSMIA", "INS", "MOTORFLUX", "RBD", "RL", "SEADL70", "SLEEP")
ContT = c("HY", "MMSE", "MOCA", "SEADL", "UPDRS_scaled", "UPDRS1_scaled", "UPDRS2_scaled", "UPDRS3_scaled", "UPDRS4_scaled")
COVs = c('FEMALE', 'YEARSEDUC', "FAMILY_HISTORY", 'AAO', "BLDfDIAG",  'DOPA', 'AGONIST')
plink_ltg = ltg_all %>% 
  mutate(FID = IID,
         IID = IID,
         FATID = 0,
         MATID = 0,
         SEX = 1+FEMALE) # In plink format, SEX 0 missing, 1 male, 2 female

models=c("TYPE;DATASET;OUTCOME;COV")

# For conditional regression, transpose observations from each individuals
trans.func = function(i){
  xi = cohort %>% filter(IID == IIDs[i]) %>% select(-IID)
  xi = as.matrix(xi)
  if(nrow(xi) >1){
    A = cumsum(rep(1, nrow(xi)))
    A1 = poly(A, degree = length(A)-1)
    transxi = t(A1) %*% xi
    transxi[abs(transxi) < 0.0000000001] = 0 # put 0
    transxi = data.frame(IID = IIDs[i], transxi)
    return(transxi)
  }else{
    return(rep(NA, length(c("OUTCOME", COVs_i_trans))+1))
  }
}

for(i in 1:length(DATASETs)){
  cat(DATASETs[i])
  # Baseline for binomial outcomes
  BL = plink_ltg %>% filter(TSTART ==0 & DATASET == DATASETs[i]) %>% 
    mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1))) # 0 missing, 1 control, 2 case
  if(DATASETs[i]=="OSLO"){ # Oslo has no appearent baseline)
    BL = plink_ltg %>% filter (COHORT == DATASETs[i]) %>% filter(-365 < TSTART & TSTART < 365) %>% 
      distinct(ID, .keep_all = T) %>% 
      mutate_at(vars(BinomT), funs(ifelse(is.na(.), 0, .+1)))
  }
      
  # PC1-5
  PC_temp = fread(paste("../dataset/", DATASETs[i], "/pca5.txt", sep = "")) 
  names(PC_temp)=c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5")
  PC = PC_temp %>% select(-FID) %>% distinct(IID, .keep_all = T) %>% mutate(IID = as.character(IID))
  
  # Covariate effective at baseline
  cov_bl = BL %>% select(FID, IID, FATID, MATID, SEX, COVs) # Cov can use 0 not as missing.
  # filter out COVARIATES with all the same values or all NAs
  temp = cov_bl[,6:ncol(cov_bl)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  COVs_bl_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  cov_bl_new = cov_bl %>% select(FID, IID, FATID, MATID, SEX, COVs_bl_i) %>% left_join(., PC, by = "IID")
    
  # Binomial outcomes
  binom_bl = BL %>% select(FID, IID, FATID, MATID, SEX, BinomT)
  # Avaial phenos for dataset_i 
  temp = binom_bl[,6:ncol(binom_bl)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  BinomTs_bl_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  if(length(BinomTs_bl_i) != 0){ # if there are such variables..
    binom_bl_new = binom_bl %>% select(FID, IID, FATID, MATID, SEX, BinomTs_bl_i)
    write.table(binom_bl_new, paste("pheno/cross/", DATASETs[i], ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    write.table(cov_bl_new, paste("pheno/cross/", DATASETs[i], "_cov.txt", sep = ""), row.names = F, quote = F, sep = "\t")
    m1=paste("cross", DATASETs[i], BinomTs_bl_i, paste(COVs_bl_i, collapse = ","), sep = ";")
  }else{BinomTs_i_cs = "";m1="None"}
  cat("    cross")

  # Continuous outcomes (need to apply different model  when longitudinal data is not available)
  cont = plink_ltg %>% filter(DATASET == DATASETs[i])
  # Covariate effective in the whole study
  cov_lt = cont %>% select(FID, IID, FATID, MATID, SEX, COVs) # Cov can use 0 not as missing.
  # filter out COVARIATES with all the same values or all NAs
  temp = cov_bl[,6:ncol(cov_bl)] %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  COVs_lt_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  cov_lt_new = cov_bl %>% select(FID, IID, FATID, MATID, SEX, COVs_lt_i) %>% left_join(., PC, by = "IID")
    
    
  # filter out OUTCOMEs with all the same values or all NAs
  temp = cont %>% select(ContT) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
  temp$NAMES = rownames(temp) # rownames were converted to the covariates.
  ContTs_i = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector # filter out NA outcomes and Var=0 outcomes
  temp2 = cont %>% filter(TSTART>0) %>% group_by(IID) %>% summarize_at(vars(ContTs_i),funs(sum(!is.na(.)))) %>% 
    summarize_at(vars(ContTs_i),funs(max)) %>% t %>% data.frame # number of observataion per ID (max)
  ContTs_i_lt=rownames(temp2)[temp2[,1] > 0] # outcome observed more than twice
  ContTs_i_bl = setdiff(ContTs_i, rownames(temp2)[temp2[,1] > 0]) # outcome observed only once.
  if(length(ContTs_i_bl) != 0){ # if there are such variables..
    cont_bl = cont %>% filter(TSTART == 0) %>% select(FID, IID, FATID, MATID, SEX, ContTs_i_bl) %>% left_join(., PC, by = "IID")
    write.table(cont_bl, paste("pheno/lnsgl/", DATASETs[i], ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    cov_bl_cont = cov_bl_new %>% rename(YEARfDIAG=BLDfDIAG) # replace BLDfDIAG -> YEARfDIAG 
    write.table(cov_bl_cont, paste("pheno/lnsgl/", DATASETs[i], "_cov.txt", sep = ""), row.names = F, quote = F, sep = "\t")
    m2=paste("lnsgl", ContTs_i_bl, paste(names(cov_bl_new)[-(1:5)], collapse = ","), sep = "_")
    cat("    lnsgl")
  }else{ContTs_i_bl = ""; m2="None";cat("    NO_lnsgl")}
  # Longitudinal set: phenotype file with covariates (for R)
  if(length(ContTs_i_lt)!=0){
    COVs_lt = gsub("BLDfDIAG", "YEARfDIAG", COVs) # replace BLDfDIAG -> YEARfDIAG
    # Filger out non-effective COVs_i_lt
    temp = cont %>% select(COVs_lt) %>% summarise_all(funs(var(.,na.rm = T))) %>% t %>% as.data.frame
    temp$NAMES = rownames(temp) # rownames were converted to the covariates.
    COVs_i_lt = temp %>% filter(V1>0) %>% select(NAMES) %>% t %>% as.vector
    cont_new = cont %>% select(FID, IID, FATID, MATID, SEX, COVs_i_lt, ContTs_i_lt) %>% left_join(., PC, by = "IID")
    write.table(cont_new, paste("pheno/lnmxs/", DATASETs[i], ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    write.table(cont_new, paste("pheno/lnmxi/", DATASETs[i], ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
    m3=paste("lnmxi", DATASETs[i], ContTs_i_lt, paste(COVs_i_lt, collapse = ","), sep = ";")
    m4=paste("lnmxs", DATASETs[i], ContTs_i_lt, paste(COVs_i_lt, collapse = ","), sep = ";")    
    cat("    lnmxi/lnmxs")
      
    # Conditinal regression model
    COVs_i_trans = COVs_i_lt[grep("^YEARfDIAG|^DOPA|^AGONIST", COVs_i_lt)]
    for (OUTCOME in ContTs_i_lt){
      cohort_temp = cont_new
      cohort_temp$OUTCOME = cohort_temp[, OUTCOME]
      cohort_temp2 = cohort_temp %>% select("IID", COVs_i_trans, "OUTCOME")
      cohort = cohort_temp2 %>% filter(!is.na(OUTCOME)) %>% arrange(IID, YEARfDIAG)
      # Transform the data to the orthogonal to the cross-sectional space
      IIDs = unique(cohort$IID)
      cohort_trans = lapply(1:length(IIDs), trans.func)
      transdata = do.call(rbind, cohort_trans)
      write.table(transdata, paste("pheno/lncns/", OUTCOME, "_", DATASETs[i], ".txt", sep=""), row.names = F, quote = F, sep = "\t")
    }
    m5=paste("lncns", DATASETs[i], ContTs_i_lt, paste(COVs_i_trans, collapse = ","), sep = ";")
    cat("    lncns")
  }else{ContTs_i_lt = ""; m3="None";m4="None";m5="None";cat("  No continous trait for longitudinal analysis")}   

  # Survival outcome
  PHENOSURV = c()
  for (j in (1:length(BinomT))){
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
      mutate(BLDfDIAG = BLDfDIAG + first(TSTART)/365.25, # Follow-up starts TSTART > 0 sometimes, Needs update
             TSTART = TSTART - first(TSTART)) %>% # New baseline set 0 (Especially, ParkWest)
      mutate(TSTOP = lead(TSTART), # The TSTART[i+1] will be the stop time
             OUTCOME = lead(OUTCOME)) %>%  # the OUTCOME is updated by the one at the stop time
      mutate(TSTOP = ifelse(is.na(TSTOP), END, TSTOP)) %>% # replace the 
      data.frame %>% arrange(IID, TSTART) %>% 
      filter(!is.na(OUTCOME)) # the basially last follow-up data will be erased (Right censored)
    
    if(nrow(cohort)!=0 & var(cohort$OUTCOME, na.rm = T) != 0){
      cohort %>% left_join(., PC, by = c("IID")) %>% 
        write.table(., paste("pheno/coxhm/", BinomT[j], "_", DATASETs[i], ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
      PHENOSURV = c(PHENOSURV, BinomT[j])
    }
  }
  m5=paste("coxhm", DATASETs[i], PHENOSURV, paste(COVs2in, collapse = ","), sep = ";") 
  cat("    coxhm \n")
  models=c(models,m1,m2,m3,m4,m5)
}
data.frame(V1=models) %>% filter(V1!="None") %>% write.table("models.txt", col.names=F, row.names = F, quote = F)
cat("all complete \n")