library(dplyr)
library(data.table)
library(metafor)

d1 = read.csv("../data/LTGall.csv")
BasicT = c('AAO', 'AD', 'FAMILY_HISTORY')
BinomT = c('HYPOSMIA', 'DEMENTIA', 'MOTORFLUX', 'DYSKINESIAS', 'DEPR', 'RL', 'HT', 'CONST', 
           'RBD', 'SLEEP', 'INS', 'HY3', "DEATH")
ContiT = c('HY', 'UPDRS1_scaled', 'UPDRS2_scaled','UPDRS3_scaled','UPDRS4_scaled','UPDRS_scaled', 
           'MMSE', 'MOCA', 'SEADL', 'PDQ39')
Covs = c('FEMALE', 'EUROPEAN', 'YEARSEDUC', 'AGEatBL', 'BLDfDIAG', 'PC1', 'PC2', 'TSTART', 'DOPA', 'AGONIST')
Cohorts = c('DATATOP', 'DIGPD', 'HBS', 'PDBP', 'PICNICS', 
            'PPMI', 'ParkWest', 'PreCEPT', 'SCOPA', 'Udall', 
            "Coriell", "ParkFit", "Predict", "Oslo") 

COEF = fread("../data/model_coeffs.csv") %>% 
  rename(cohort=V2, outcome=V3, predictor=V4, variable = V1, mean=V5, se=V6, pv=V8, NinSet=V9) %>% 
  mutate_at(vars(c('mean', 'se', 'pv')), as.numeric) %>% 
  filter(!is.na(mean))

COEF %>%  filter(cohort == "DATATOP" & outcome =="AAO" & variable =="PC1") %>% with(hist(pv))
COEF %>%  filter(cohort == "DATATOP" & outcome =="AAO" & variable =="PC2") %>% with(hist(pv))
COEF %>%  filter(cohort == "PPMI" & outcome =="DEMENTIA" & model == "coxph" & variable =="YEARSEDUC") %>% with(hist(pv))

COEF2 = COEF %>%
  filter(!variable == "(Intercept)") %>% 
  filter(!variable == "PREDICTOR") %>% 
  group_by(cohort, outcome, model, variable) %>% 
  summarize_at(vars("mean", "se", "pv"), median) %>% 
  data.frame

models = unique(COEF2$model)
outcomes = unique(COEF2$outcome)
variables = unique(COEF2$variable)
x = 0
meta_result = data.frame(matrix(NA, nrow=1000, ncol=10))
names(meta_result) = c('model', 'outcome', 'variable', 'Method', 'NoCohort', 'meta_beta', 'meta_se',
                       'pvalue', 'Homogeneity', 'I2')
for(i in 1:length(models)){
  for (j in 1:length(outcomes)){
    for (k in 1:length(variables)){
      temp = COEF2 %>% 
        filter(model == models[i] & outcome == outcomes[j] & variable == variables[k])
      print(nrow(temp))
      if(nrow(temp)==0){next}
      res_fe = try(rma(yi = mean, sei = se,  data=temp, method = "FE"), silent = T)
      if(class(res_fe)[1]=="try-error"){
        meta_result[x,] = 
          c(models[i], outcomes[j], variables[k], "NoConverge", rep(NA,6))
        x = x +1
        next
      }
      g = function(){
        res = try(rma(yi = mean, sei = se,  data=temp, method = "REML"),silent = T)
        # Some can't converge by REML so use alternative method
        if(class(res)[1]=='try-error'){
          res = try(rma(yi = mean, sei = se,  data=temp, method = "ML"), silent = T)
          if(class(res)[1]=='try-error'){
            res = data.frame(method = "NonConML", QEp="NA", I2 = "NA")
          }
          method = 'ML'
        }
        res
      }
      res_me = g()
      meta_result[x,] = 
        c(models[i], outcomes[j], variables[k], res_me$method, nrow(temp),
          res_fe$beta, res_fe$se, res_fe$pval, res_me$QEp, res_me$I2)
      x=x+1
    }
  }
}
meta_result = meta_result %>% filter(!is.na(NoCohort))
meta_result[,5:9]=lapply(meta_result[,5:9], as.numeric)
temp = meta_result %>% 
  filter(pvalue<0.05) %>% 
  filter(I2<80) %>% 
  arrange(model, outcome, pvalue)

temp2 = meta_result %>% 
  filter(pvalue>0.8) %>% 
  filter(I2<10) %>% 
  arrange(model, outcome, pvalue)
