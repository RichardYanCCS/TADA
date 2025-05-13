rm(list=ls())
options(mc.cores = parallel::detectCores()) 

library(dplyr)
library(tidyr)
library(readr)
library(haven)

# load in subject-level dataset
# adtte <- haven::read_sas("./sourceDataIPD/sasDataSet/adtte.sas7bdat") 
# saveRDS(adtte, file = "./sourceDataIPD/sasDataSet/adtte.RDS")

# data loading takes a while, use the saved RDS
adsl <- haven::read_sas("./sourceDataIPD/sasDataSet/adsl.sas7bdat")
adtte <- readRDS("./sourceDataIPD/sasDataSet/adtte.RDS")

# OS and CNSR (censor indicator)
subsetOS <- adtte %>% 
  filter(PARAM == "Overall Survival Time (months)") %>% 
  select(USUBJID, PARAM, AVAL, EVNTDESC, CNSR)

# extend IPCW covariates: AGE, SEXN, ECOGBL, RACEN, SMKBN
IPD <- bind_cols(select(adsl, AGE, SEXN, RACEN, SMKBN, ECOGBL), 
                 select(subsetOS, AVAL, EVNTDESC, CNSR)) %>%
  mutate(RACEN = replace(RACEN, RACEN == 4, 6)) %>% 
  mutate(SEX = factor(SEXN,  levels = c(1, 2), labels = c("Male", "Female")), 
         RACE = factor(RACEN, levels = c(5, 3, 6, 2), labels = c("White", "Black or African American", "Other", "Asian")),
         ECOGBL = factor(ECOGBL, levels = c(0, 1, 2, 3), labels = c("0", "1", "2", "3")),
         SMOKE = factor(SMKBN, levels = c(1, 2, 3), labels = c("Light EX-Smoker", "Non-Smoker", "Smoker")),
         OS = AVAL, 
         DEATH = 1-CNSR, 
         CENS = CNSR,
         EVENT_REASON = EVNTDESC) %>%
  select(-c(SEXN, RACEN, SMKBN, AVAL, EVNTDESC, CNSR)) %>%
  drop_na()


# IPD data saving
output_IPD_path <- "./processed_source_IPD.csv"
write_csv(IPD, output_IPD_path)

# AgD data generation, based on Table 1
AgD <- data.frame(2056, 72, 0.651) 
colnames(AgD) <- c("N", "AGE_MEAN", "SEX_PROP")

output_AgD_path <- "./processed_target_AgD.csv"  
write_csv(AgD, output_AgD_path)

