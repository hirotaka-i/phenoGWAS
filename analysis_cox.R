# Cox model for GWAS
library(data.table)
library(dplyr)
library(zoo)
library(survival)
library(parallel)


# Create temp data
temp =data.frame(ID = paste(1:5000, 1:5000, sep="_"))
temp[, 2:10001] = rnorm(5000*10000, 0, 1)
write.table(temp, "temp.txt", row.names = F, quote = F, sep = "\t")




# Read data
FILES = list.files("outputs/surv/",full.names = T)
EVAL = data.table(NAME = as.character(FILES))
EVAL$STUDY = lapply(strsplit(FILES, "\\."), "[", 2) %>% unlist
EVAL$OUTCOME = lapply(strsplit(FILES, "\\."), "[",1) %>% unlist %>% substring(., 15)
EVAL = EVAL %>% arrange(STUDY, OUTCOME)
DATASETs= unique(EVAL$STUDY)

for (i in 1:length(DATASETs)){ # iterate cohort
  SET = EVAL %>% filter(STUDY == DATASETs[i])
  
  # read SNP data
  t = Sys.time()
  SNPset = fread("temp.txt")
  SNPs = names(SNPset)[2:ncol(SNPset)]
  
  for (j in 1:nrow(SET)){ # iterate OUTCOMEs
    cohort=fread(SET$NAME[j]) %>% arrange(ID, TSTART)
    cohort_snp = left_join(cohort, SNPset, by = "ID")
    cohort_snp$SurvObj1 = with(cohort_snp, Surv(TSTART, TSTOP, OUTCOME == 1))

    surv.func2 = function(x){
      # Models
      MODEL5 = paste("SurvObj1 ~", paste(c(SNPs[x], names(cohort)[4:(length(names(cohort))-8)], paste("PC",1:5, sep="")), collapse = "+") )
      testCox5 = try(coxph(eval(parse(text = MODEL5)), data = cohort_snp),silent = T)
      if(class(testCox)[1]=='try-error'){next}
      RES = summary(testCox5)$coefficients[1,]
      EVENT_OBS = paste(testCox5$nevent, testCox5$n, sep="_")
      sumstat <- c(SNPs[x], EVENT_OBS, as.numeric(RES[4]), RES[1], RES[3], RES[5])
    }
    temp = mclapply(1:length(SNPs), surv.func2)
    temp2 = do.call(rbind, temp)
    attributes(temp2)$dimnames[[2]]=c("SNP", "EVENT_OBS", "TEST", "Beta", "SE", "Pvalue")
    write.table(temp2, paste(SET$OUTCOME[j], DATASETs[i], "chrtemp.txt", sep="_"), row.names = F, quote = F, sep = "\t")
  }
  
  Sys.time()-t
}


rm(list=ls())
