# Cox model for GWAS

rm(list=ls())
Sys.setlocale('LC_ALL', locale='C')
# setwd("/Users/iwakih2/Downloads/META_new")


library(data.table)
library(dplyr)
library(tableone)
library(zoo)
library(survival)


##### FACOTORS ######################################################################################
Covs = c('FEMALE', 'YEARSEDUC', 'AGEatBL', 'BLDfDIAG','TSTART', 'DOPA', 'AGONIST')
BasicT = c('AAO', 'FAMILY_HISTORY')
BinomT = c('HYPOSMIA', 'DEMENTIA', 'MOTORFLUX', 'DYSKINESIAS', 'DEPR', 'RL', 'CONST', 'RBD', 'SLEEP', 'INS', 'HY3')
# DEATH often happened after the observation period so should be analyzed separately. 
predictors <- c("grsZ","NeuroX_dbSNP_rs114138760_C","exm106217_C","exm106220_A","NeuroX_rs2230288_A",
                "NeuroX_rs823118_C","NeuroX_rs10797576_T","NeuroX_rs6430538_T","NeuroX_rs1955337_T",
                "NeuroX_rs12637471_A","NeuroX_dbSNP_rs34884217_replciate_1_G","NeuroX_rs34311866_G",
                "NeuroX_rs11724635_C","exm.rs6812193_T","NeuroX_rs356181_T","NeuroX_rs3910105_C",
                "exm535099_T","NeuroX_dbSNP_rs115462410_T","NeuroX_rs199347_C","NeuroX_rs591323_A",
                "NeuroX_dbSNP_rs118117788_T","NeuroX_rs329648_T","NeuroX_rs76904798_T","exm994671_A",
                "NeuroX_rs11060180_G","NeuroX_rs11158026_T","NeuroX_rs2414739_G",
                "NeuroX_dbSNP_rs14235_replciate_1_A","exm.rs11868035_A","NeuroX_rs17649553_T",
                "NeuroX_rs12456492_G","NeuroX_rs55785911_A")
PCs = c("PC1", "PC2")
######################################################################################################
ltg_all=fread('../data/LTGall.csv') %>% data.frame %>%
  filter(STUDY_NAME != "Predict") %>% 
  rename(COHORT = STUDY_NAME) %>% 
  mutate(COHORT= case_when(
    COHORT == "Coriell" ~ "NET-PD_LS1",
    COHORT == "SCOPA" ~ "PROPARK",
    COHORT == "PreCEPT" ~ "PRECEPT",
    COHORT == "ParkFit" ~ "PARKFIT",
    COHORT == "ParkWest" ~ "PARKWEST",
    COHORT == "Udall" ~ "UDALL",
    TRUE ~ COHORT
  )) %>% 
  select(c("COHORT", "ID", Covs, BasicT, BinomT, PCs, predictors))

Oslo = HYUPDRSg %>% 
  mutate(COHORT = "Oslo") %>% 
  select(names(ltg_all))

data_all = rbind(ltg_all, Oslo) %>% 
  arrange(COHORT)

N_VARS = 2 + length(c(Covs, BasicT, BinomT, PCs))

#######################################################################
# models specifications
COHORTS = unique(data_all$COHORT)
COHORTS
MODELS = data.frame(
  COHORT = COHORTS,
  MODEL = c(
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE                              + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE             + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR             + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG"
  )
)
# Another way is to replace missing covariates with the same value like "MISS"
GRID = expand.grid(1:length(BinomT), 1:length(predictors))
cox.stats.func <- function(x){
  STUDY = COHORTS[i]
  OUTCOME = BinomT[GRID[x,1]]
  PREDICTOR = predictors[GRID[x,2]]
  # print(paste(i, x, STUDY, OUTCOME, PREDICTOR))
  ltg = subset(data_all, COHORT == STUDY)
  ltg$OUTCOME <- ltg[,OUTCOME]
  ltg$PREDICTOR = ltg[,PREDICTOR]
  ltg = ltg %>% select(-one_of(predictors))
  
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
  if(nrow(cohort)==0){
    sumstat <- cbind(STUDY, OUTCOME, PREDICTOR, "NoData", "NoData", "NoData", "NoData", "NoData")
  }else if(var(cohort$OUTCOME, na.rm = T) == 0 | var(cohort$PREDICTOR, na.rm = T) == 0){
    sumstat <- cbind(STUDY, OUTCOME, PREDICTOR, "NoVar", "NoVar", "NoVar", "NoVar", "NoVar")
  }else{
    cohort$SurvObj1 = with(cohort, Surv(TSTART, TSTOP, OUTCOME == 1))
    testCox = try(coxph(eval(parse(text = paste(MODELS[i,2]))), data = cohort),silent = T)
    # Some can't converge by REML so use alternative method
    if(class(testCox)[1]=='try-error'){
      sumstat <- cbind(STUDY, OUTCOME, PREDICTOR, "OverFlow", "OverFlow", "OverFlow", "OverFlow", "OverFlow")
    }else{
      temp = summary(testCox)$coefficients
      betas <- temp[,1]
      se <- temp[,3]
      zval <- temp[,4]
      pval <- temp[,5]
      NinSet = paste('E', testCox$nevent, '/N', testCox$n, sep='')
      sumstat <- cbind(STUDY, OUTCOME, PREDICTOR, betas, se, zval, pval, NinSet)
    }
  }
  sumstat
}

temp = mclapply(1:nrow(GRID), cox.stats.func, mc.cores = 2)
temp2 = do.call(rbind, temp)

# t <- proc.time()
# invisible(mclapply(1:nrow(GRID), cox.stats.func, mc.cores = 2))
# duration3 <- (proc.time() - t)[3]
# t <- proc.time()
# invisible(lapply(1:nrow(GRID), cox.stats.func))
# duration2 <- (proc.time() - t)[3]
