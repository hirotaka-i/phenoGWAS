# Cox model for GWAS

library(data.table)
library(dplyr)
library(zoo)
library(survival)
library(parallel)
# Read data
ltg_all=fread('outputs/PhenoFile.csv') %>% data.frame
names(ltg_all)

# Factors
PCs = c("PC1", "PC2")
COVs = c("ID", 'FEMALE', 'YEARSEDUC', "FAMILY_HISTORY", 'AAO', "BLDfDIAG",  'DOPA', 'AGONIST', "TSTART")
OUTCOMEs = c('HYPOSMIA', 'DEMENTIA', 'MOTORFLUX', 'DYSKINESIAS', 'DEPR', 'RL', 'CONST', 'RBD', 'SLEEP', 'INS', 'HY3')
COHORTs = unique(ltg_all$COHORT)
N_VARS = 2 + length(c(Covs, BasicT, BinomT, PCs))
PREDICTORs = paste("X", 1:100, sep = "")
#######################################################################
# models specifications
COHORTs
MODELs = data.frame(
  COHORT = COHORTs,
  MODEL = c(
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE             + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2                  + FEMALE                              + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG",
    "SurvObj1 ~ PREDICTOR + PC1 + PC2 + DOPA + AGONIST + FEMALE + YEARSEDUC + FAMILY_HISTORY + AAO + BLDfDIAG"
    )
)
MODELs

# Another way is to replace missing covariates with the same value like "MISS"
ncol(ltg_all)

# For test
ltg_all = ltg_all %>% mutate(PC1 = rnorm(nrow(ltg_all), 0, 1), 
                             PC2 = rnorm(nrow(ltg_all), 0, 1))
ltg_all[,42:141] = rnorm(nrow(ltg_all)*100, 0, 1)

PREDICTORs = paste("V", 42:141, sep = "")
# Survival Analaysis
GRID = expand.grid(1:length(OUTCOMEs), 1:length(PREDICTORs))
cox.stats.func <- function(x){
  STUDY = COHORTs[i]
  OUTCOME = OUTCOMEs[GRID[x,1]]
  PREDICTOR = PREDICTORs[GRID[x,2]]
  # print(paste(i, x, STUDY, OUTCOME, PREDICTOR))
  ltg = ltg_all %>% filter(COHORT == STUDY)
  ltg$OUTCOME = ltg[,OUTCOME]
  ltg$PREDICTOR = ltg[,PREDICTOR]
  ltg = ltg %>% select(OUTCOME, PREDICTOR, COVs, PCs)
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
    testCox = try(coxph(eval(parse(text = paste(MODELs[i,2]))), data = cohort),silent = T)
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
  sumstat[1,]
}

t <- proc.time()
i = 9
temp = mclapply(1:nrow(GRID), cox.stats.func)
temp2 = do.call(rbind, temp)
proc.time() - t


rm(list=ls())
