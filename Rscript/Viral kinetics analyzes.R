##################################################################################
# Author: Ilse Westerhof
# Last updated: 19 AUG 2024 
# Program: MACOS 14.3.1, R version 4.4.1, Rstudio version 2024.04.2+764
# https://github.com/IWesterhof/Antibody-profiles-and-COVID-19

##################################################################################

########## Prep session ########## 
rm(list = ls())
Sys.setlocale(locale="English")

## Load  R packages
pacman::p_load(readr, tidyverse, magrittr,  dplyr, devtools, reshape2,  ggplot2, gtsummary, arsenal, cowplot)

########## Prep session ########## 
rm(list = ls())
Sys.setlocale(locale="English")

## Work directories
base.dir    <- "~/Networkshares/datamanagement/Research/ID/VERDI/E_ResearchData/0_DataPreparation/Rproject Antibodies and disease burden"
data.dir    <- file.path(base.dir, "01_Data")
code.dir    <- file.path(base.dir, "02_Code")
output.dir  <- file.path(base.dir, "03_Output")

getwd()
setwd("~/Networkshares/datamanagement/Research/ID/VERDI/E_ResearchData/0_DataPreparation/Rproject Antibodies and disease burden")

### Load data
load("01_Data/VERDI.RData"); load("01_Data/VERDIpp.RData"); load("01_Data/AllParticipants.RData")
Samples_VERDI <- read_csv("~/Networkshares/datamanagement/Research/ID/VERDI/E_ResearchData/3_Master/Samples_VERDI.csv", col_types = cols(...1 = col_skip()))
Samples_VERDI_Coprimary <- Samples_VERDI %>% select(ref_id, sample_type, collection_date, result) %>%
  filter(!sample_type=="dbs") %>% arrange(ref_id, collection_date, desc(result)) %>% distinct(ref_id, sample_type, .keep_all=T) %>%
  ungroup() %>% arrange(ref_id, desc(result)) %>% distinct(ref_id, .keep_all=T) %>% rename(coprimary = result) %>% select(ref_id, coprimary)
hhm_refids <- read_csv("01_Data/hhm_refids.csv") %>% mutate(IncludedAntibodies="Yes")
hhm_refids = merge(hhm_refids, Samples_VERDI_Coprimary, by="ref_id", all=T)  %>% filter(IncludedAntibodies=="Yes")
AgeIndex <- AllParticipants %>% filter(Index=="index") %>% select(household_id, Age_group) %>% distinct() %>% 
  filter(!household_id=="nl082") %>% rename(Age_group_index = Age_group) # removed because not age of index is unknown
PriorInfection <- AllParticipants %>% filter(Index=="household member") %>% select(ref_id, PriorInfection) %>% distinct() 

#### subset positive household members #### 
VERDIpp <- merge(VERDIpp, hhm_refids, by="ref_id", all=T) %>% filter(IncludedAntibodies=="Yes") %>% select(-contains("_60000"))
VERDI <- merge(VERDI, hhm_refids, by="ref_id", all=T) %>% filter(IncludedAntibodies=="Yes")
AllParticipants_Infection  = AllParticipants %>% filter(Index == "household member") %>%   filter(!(is.na(Age) | is.na(`NL63-S1-T_start`) | is.na(SEC_TRANS))) 
AllParticipants <- merge(AllParticipants, hhm_refids, by="ref_id", all=T) %>% filter(IncludedAntibodies=="Yes")
Samples_VERDI <- merge(Samples_VERDI, hhm_refids, by="ref_id", all=T) %>% filter(IncludedAntibodies=="Yes")

classify_based_on_median <- function(column) {
    median_value <- median(column, na.rm = TRUE)  # Calculate median, handling NAs if present
    ifelse(column < median_value, "low", "high")}
antibodies <- VERDIpp %>%
  mutate(across(c(`SARS2-ecto-T_start`, `SARS-CoV-2 WT S1_start`, `SARS-CoV-2 WT np_start`, `NL63-S1-T_start`, `NL63-ecto_start`, 
                  `229E-S1_start`, `229E-NP_start`, `HKU1-S1-mFc_start`, `HKU1-NP_start`, `OC43-ecto_start`),    ~ classify_based_on_median(.), .names = "{sub('_start$', '_median', .col)}"))  %>% 
  select(ref_id, contains("median"), Age) %>% select(-`SARS-CoV-2 WT RBD_median`) %>% 
  rename(Influenza_median = `Influenza A H1N1 (A/California/04/2009)_median`)
  
sympt = VERDIpp %>% select(ref_id, Diseaseseverity, onsetsympt, endsympt, StartEpi, EndEpi)
Viralkinetics = Samples_VERDI %>% select(ref_id, household_id) %>% distinct()

#### Determine viral kinetic parameters #### 
Samples_VERDI <- merge(Samples_VERDI, antibodies, by="ref_id") %>% filter(sample_type=="nts" | sample_type=="saliva") %>% mutate(result_ct = ifelse(is.na(result_ct) | result_ct==0, 40, result_ct)) %>%
  select(ref_id, samplenumber, sample_type, collection_date, result_ct)
peak = Samples_VERDI %>% filter(result_ct<40) %>% group_by(ref_id, sample_type) %>% mutate(n_collectedsamples = length(sample_type)) %>% 
    arrange(result_ct)   %>%  select(collection_date, n_collectedsamples, ref_id, sample_type, result_ct) %>% ungroup() %>% distinct(ref_id, sample_type, .keep_all=T) %>% rename(Peak=collection_date, Peak_ct=result_ct) ; peak
Firstbelow40 = Samples_VERDI %>% filter(result_ct<40) %>% arrange(ref_id, collection_date) %>% group_by(ref_id, sample_type) %>% distinct(ref_id, sample_type, .keep_all=T) %>% select(collection_date) %>% rename(Firstbelow40=collection_date);Firstbelow40
Lastbelow40 = Samples_VERDI %>% filter(result_ct<40) %>% arrange(ref_id, desc(collection_date)) %>% group_by(ref_id, sample_type) %>% distinct(ref_id, sample_type, .keep_all=T) %>% select(collection_date) %>% rename(Lastbelow40=collection_date);Firstbelow40
Firstbelow30 = Samples_VERDI %>% filter(result_ct<30) %>% arrange(ref_id, collection_date) %>% group_by(ref_id, sample_type) %>% distinct(ref_id, sample_type, .keep_all=T) %>% select(collection_date) %>% rename(Firstbelow30=collection_date);Firstbelow40
Lastbelow30 = Samples_VERDI %>% filter(result_ct<30) %>% arrange(ref_id, desc(collection_date)) %>% group_by(ref_id, sample_type) %>% distinct(ref_id, sample_type, .keep_all=T) %>% select(collection_date) %>% rename(Lastbelow30=collection_date);Firstbelow40
Viralkinetics = merge(Viralkinetics, peak, by=c("ref_id"), all=T)
Viralkinetics = merge(Viralkinetics, Firstbelow40, by=c("ref_id", "sample_type"), all=T)
Viralkinetics = merge(Viralkinetics, Lastbelow40, by=c("ref_id", "sample_type"), all=T)
Viralkinetics = merge(Viralkinetics, Firstbelow30, by=c("ref_id", "sample_type"), all=T)
Viralkinetics = merge(Viralkinetics, Lastbelow30, by=c("ref_id", "sample_type"), all=T)
Viralkinetics = Viralkinetics %>% mutate(durationbelow40 = Lastbelow40-Firstbelow40+1,
                                         durationbelow30 = Lastbelow30-Firstbelow30+1)

Viralkinetics = merge(Viralkinetics, antibodies, by="ref_id", all=T)
Viralkinetics = merge(sympt, Viralkinetics, by="ref_id", all=T) 
Viralkinetics = merge(Viralkinetics, AgeIndex, by="household_id", all=T) %>%
  mutate(Cum_alpha = rowSums(select(., contains(c("229E", "NL63"))) == 'high', na.rm=T),
         Cum_beta = rowSums(select(., contains(c("SARS2", "SARS-CoV-2", "HKU1", "OC43"))) == 'high', na.rm=T),
         Cum_SARS2 = rowSums(select(., contains(c("SARS2", "SARS-CoV-2"))) == 'high', na.rm=T),
         Cum_NL63 = rowSums(select(., contains("NL63")) == 'high', na.rm=T),
         Cum_229E = rowSums(select(., contains("229E")) == 'high', na.rm=T),
         Cum_HKU1 = rowSums(select(., contains("HKU1")) == 'high', na.rm=T),
         Cum_OC43 = rowSums(select(., contains("OC43")) == 'high', na.rm=T))
Samples_VERDI = merge(Samples_VERDI, Viralkinetics, by=c("ref_id","sample_type"), all=T) 
Samples_VERDI = Samples_VERDI %>% mutate(daytopeak = collection_date-Peak) %>%
  mutate(Cum_alpha = rowSums(select(., contains(c("229E", "NL63"))) == 'high', na.rm=T),
         Cum_beta = rowSums(select(., contains(c("SARS2", "SARS-CoV-2", "HKU1", "OC43"))) == 'high', na.rm=T),
         Cum_SARS2 = rowSums(select(., contains(c("SARS2", "SARS-CoV-2"))) == 'high', na.rm=T),
         Cum_NL63 = rowSums(select(., contains("NL63")) == 'high', na.rm=T),
         Cum_229E = rowSums(select(., contains("229E")) == 'high', na.rm=T),
         Cum_HKU1 = rowSums(select(., contains("HKU1")) == 'high', na.rm=T),
         Cum_OC43 = rowSums(select(., contains("OC43")) == 'high', na.rm=T))

#### Table 1. Characteristics of SARS-CoV-2 exposed household members. ####
table1data = merge(AllParticipants_Infection, AgeIndex, by="household_id", all=T) %>%  
  select(ref_id, Age, Age_group_index, Vaccinated, contains("_start"), SEC_TRANS) %>%
  filter(!is.na(`SARS-S1_start`)) 
table1data = merge(table1data, PriorInfection, by="ref_id", all=T) 
table1data %>% select(SEC_TRANS, Age, Vaccinated, PriorInfection, Age_group_index, contains("_start")) %>%
  select(-starts_with("B_"), -contains(c("MERS", "Influenza", "RBD", "SARS-ecto", "SARS-S1", "SARS-NP"))) %>% 
  tbl_summary(by = SEC_TRANS) %>%  add_overall()

#### Table 2. SARS-CoV-2 infection risk by baseline antibody status for different coronaviruses; aORs for getting infected per one level increase in cumulative antibody score. ####
Coprimary = AllParticipants %>% select(ref_id, coprimary)
Infection = merge(AllParticipants_Infection, AgeIndex, by="household_id", all=T) %>%  
  select(ref_id, Age, Age_group_index, contains("_start"), SEC_TRANS) %>%
  filter(!is.na(`SARS-S1_start`)) %>%
  mutate(across(c(`SARS2-ecto-T_start`, `SARS-CoV-2 WT S1_start`, `SARS-CoV-2 WT np_start`, `NL63-S1-T_start`, `NL63-ecto_start`, 
                  `229E-S1_start`, `229E-NP_start`, `HKU1-S1-mFc_start`, `HKU1-NP_start`, `OC43-ecto_start`),  ~ classify_based_on_median(.), .names = "{sub('_start$', '_median', .col)}")) %>%
  mutate(Cum_alpha = rowSums(select(., contains(c("229E", "NL63"))) == 'high', na.rm=T),
         Cum_beta = rowSums(select(., contains(c("SARS2", "SARS-CoV-2", "HKU1", "OC43"))) == 'high', na.rm=T),
         Cum_SARS2 = rowSums(select(., contains(c("SARS2", "SARS-CoV-2"))) == 'high', na.rm=T),
         Cum_NL63 = rowSums(select(., contains("NL63")) == 'high', na.rm=T),
         Cum_229E = rowSums(select(., contains("229E")) == 'high', na.rm=T),
         Cum_HKU1 = rowSums(select(., contains("HKU1")) == 'high', na.rm=T),
         Cum_OC43 = rowSums(select(., contains("OC43")) == 'high', na.rm=T),
         SEC_TRANS2 = ifelse(SEC_TRANS == "Secondary transmission", 1, 0)) %>% 
  filter(!is.na(Age_group_index))
Infection = merge(Infection, PriorInfection, by="ref_id", all=T) 
Coprimary = merge(Coprimary, Infection, by="ref_id", all=T) %>% filter(!coprimary=="Positief")

Infection %>% filter(!is.na(Age)) %>% select(contains("Cum_"), "SEC_TRANS") %>% 
  mutate(SEC_TRANS = factor(SEC_TRANS, levels = c("Secondary transmission", "No secondary transmission"))) %>% 
  tbl_summary(by = SEC_TRANS, missing = "no",
   type = all_of(contains("Cum_")) ~ "continuous2", 
   statistic = all_continuous() ~ c("{median} ({p25}, {p75})"), 
   digits = all_continuous() ~ 0) %>%  add_overall()

tbl_regression(glm(SEC_TRANS2~Cum_alpha + Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~Cum_beta + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~SARS2_Ecto + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~SARS2_S1 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~SARS2_NP + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~ Cum_NL63 + Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~ Cum_229E + Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~ Cum_HKU1 + Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 
tbl_regression(glm(SEC_TRANS2~ Cum_OC43 + Cum_SARS2 + Age + Age_group_index, data=Infection, family = binomial(link = logit)),exponentiate = TRUE) 

#### Table 3. Symptomatology of SARS-CoV-2 infection by baseline antibody status for different coronaviruses ####
sub = merge(AllParticipants, AgeIndex, by="household_id", all=T) %>% 
  filter(!(is.na(Age_group_index) | is.na(IncludedAntibodies))) %>%
  select(ref_id, Age,Age_group_index, contains("_start"), Diseaseseverity, DaysWSympt, cumSeverityscore) %>%
  mutate(across(c(`SARS2-ecto-T_start`, `SARS-CoV-2 WT S1_start`, `SARS-CoV-2 WT np_start`, `NL63-S1-T_start`, `NL63-ecto_start`, 
                  `229E-S1_start`, `229E-NP_start`, `HKU1-S1-mFc_start`, `HKU1-NP_start`, `OC43-ecto_start`),  ~ classify_based_on_median(.), .names = "{sub('_start$', '_median', .col)}")) %>%
  mutate(Cum_alpha = rowSums(select(., contains(c("229E", "NL63"))) == 'high', na.rm=T),
         Cum_beta = rowSums(select(., contains(c("SARS2", "SARS-CoV-2", "HKU1", "OC43"))) == 'high', na.rm=T),
         Cum_SARS2 = rowSums(select(., contains(c("SARS2", "SARS-CoV-2"))) == 'high', na.rm=T),
         Cum_NL63 = rowSums(select(., contains("NL63")) == 'high', na.rm=T),
         Cum_229E = rowSums(select(., contains("229E")) == 'high', na.rm=T),
         Cum_HKU1 = rowSums(select(., contains("HKU1")) == 'high', na.rm=T),
         Cum_OC43 = rowSums(select(., contains("OC43")) == 'high', na.rm=T)) %>%
  select(ref_id, Age,Age_group_index, contains("Cum"), Diseaseseverity, DaysWSympt, cumSeverityscore, contains("63")) 

model = nnet::multinom(Diseaseseverity ~ Cum_SARS2 + Age + Age_group_index, data = sub)
model2 = nnet::multinom(Diseaseseverity ~ Age + Age_group_index, data = sub)
tbl_regression(model, exponentiate=T); lr_test <- lmtest::lrtest(model, model2); cat("\nOverall p-value:", round(lr_test$Pr[2], digits = 2))

tbl_regression(lm(cumSeverityscore~Cum_beta + Age + Age_group_index, data=sub))
tbl_regression(lm(DaysWSympt~Cum_beta + Age + Age_group_index, data=sub))


#### Table 4. CT-value trajectories in SARS-CoV-2 infected subjects by baseline antibody status for different coronaviruses; adjusted mean difference per one level increase in cumulative antibody score.  ####
nts <- Viralkinetics %>% filter(sample_type=="nts") #
nts%>% select(contains(c("Peak_ct", "durationbelow")), sample_type) %>%
  tbl_summary(by = sample_type, missing = "no", statistic = all_continuous() ~ c("{mean} ({sd})")) 
tbl_regression(glm(as.numeric(durationbelow40)           ~ Cum_beta + Age + Age_group_index, data=nts)) 
tbl_regression(glm(as.numeric(durationbelow30)           ~ Cum_beta + Age + Age_group_index, data=nts)) 
tbl_regression(glm(as.numeric(Peak_ct)                   ~ Cum_beta + Age + Age_group_index, data=nts)) 


#### Figure 1 CT-value trajectories for SARS-CoV-2 infected individuals. The observed viral load peak is set at day 0 for each subjects   #### 
Agectvalue= Samples_VERDI %>% filter(!is.na(`SARS2-ecto-T_median`)& !is.na(sample_type)) %>% 
  mutate(AgeGroup = ifelse(Age<12, "Child",  ifelse(Age>=12 & Age<18, "Adolescent",  ifelse(Age>=18, "Adult", NA)))) %>%
  mutate(AgeGroup = factor(AgeGroup, levels=c("Child", "Adolescent", "Adult"))) %>%
  ggplot(aes(daytopeak, result_ct, colour=AgeGroup)) +geom_vline(xintercept=0) +  geom_point(alpha=0.4) + geom_smooth() + 
  scale_x_continuous(limits = c(-15, 20.1), expand = c(0,0))  + 
  scale_y_reverse(limits = c(40.2,10), expand = c(0,0))  + scale_colour_brewer(palette = "Set1") + labs(colour="Age", x= "Day to peak", y="Ct value") +
  facet_grid(.~sample_type)  + scale_colour_brewer(palette = "Set1") +theme_bw()+theme(legend.position = "bottom")
Severity = Samples_VERDI %>% filter(!is.na(`SARS2-ecto-T_median`)& !is.na(sample_type)) %>% 
  mutate(Diseaseseverity = ifelse(Diseaseseverity=="ARI", "Symptomatic",  ifelse(Diseaseseverity=="Mild symptoms", "Pauci-symptomatic",  ifelse(Diseaseseverity=="Asymptomatic", "Asymptomatic", Diseaseseverity)))) %>%
  mutate(Diseaseseverity = factor(Diseaseseverity, levels=c("Symptomatic", "Pauci-symptomatic", "Asymptomatic"))) %>%
  ggplot(aes(daytopeak, result_ct, colour=Diseaseseverity)) +geom_vline(xintercept=0) +  geom_point(alpha=0.4) + geom_smooth() + 
  scale_x_continuous(limits = c(-15, 20.1), expand = c(0,0))  + 
  scale_y_reverse(limits = c(40.2,10), expand = c(0,0))  + scale_colour_brewer(palette = "Set1") + labs(colour="Disease severity", x= "Day to peak", y="Ct value") +
  facet_grid(.~sample_type)  + scale_colour_brewer(palette = "Set1") +
plot_grid(Agectvalue, Severity, align = "hv", labels = "AUTO", ncol = 1)

pdf(file = "/Users/iwester4/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Antibody profiles/Manuscript/Figure 1 CT-value trajectories for SARS-CoV-2 infected individuals.pdf", width = 10.5, height = 8.5)
plot_grid(Agectvalue, Severity, align = "hv", labels = "AUTO", ncol = 1)
dev.off()

#### Figure 2 CT-value trajectories for SARS-CoV-2 infected individuals by baseline antibody status for different coronaviruses (SARS-CoV-2, NL63, 229E, HKU1, OC43).   #### 
pdf(file = "/Users/iwester4/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Antibody profiles/Manuscript/Figure 2 CT-value trajectories for SARS-CoV-2 infected individuals by baseline antibody status for different coronaviruses.pdf", width = 10.5, height = 5.5)
Samples_VERDI %>% filter(!is.na(`SARS2-ecto-T_median`) & sample_type=="nts") %>% 
  dplyr::select(-contains(c("B_", "MERS", "Influenza", "SARS-ecto_start", "SARS-S1_start", "SARS-CoV-2 WT RBD_start", "SARS-NP_start"))) %>% 
  pivot_longer(contains("_median"), names_to = "Virus", values_to = "value") %>%  filter(!is.na(value)) %>%
  mutate(target = ifelse(grepl("S1", Virus), "S1", ifelse(grepl("NP|np", Virus), "NP","Ecto"))) %>% 
  mutate(Virus = str_replace(Virus, "_start", "") ,Virus = str_extract(Virus, "^[^-]+")) %>% mutate(Virus = ifelse(Virus=="SARS", "SARS2", Virus)) %>% 
  mutate(Virus = factor(Virus, levels=c("SARS2", "NL63", "229E", "HKU1", "OC43"))) %>%
  ggplot(aes(daytopeak, result_ct, colour=value)) +geom_vline(xintercept=0) +  geom_point(alpha=0.4) + geom_smooth() + 
  scale_x_continuous(limits = c(-15, 20.1), expand = c(0,0))  + 
  scale_y_reverse(limits = c(40.2,10.1), expand = c(0,0))  + scale_colour_brewer(palette = "Set1") + labs(colour="Antibody titer", x= "Day to peak", y="Ct value") +
  facet_grid(target~Virus)  + scale_colour_brewer(palette = "Set1") +theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), strip.text = element_text(face = "bold")) 
dev.off()

#### Supplementary Fig. 1 Correlation of cumulative high antibody levels for different coronaviruses measured at baseline. #### 
library(PerformanceAnalytics)
pdf(file = "/Users/iwester4/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Antibody profiles/Manuscript/Supplementary Figure 1 Correlation of cumulative high antibody levels for different coronaviruses measured at baseline.pdf", width =7, height = 6.5)
sub <- Viralkinetics %>% dplyr::select(contains("Cum")) 
colnames(sub) <- gsub("^Cum_", "", colnames(sub))
sub <- sub %>% dplyr::select(-alpha, -beta) 
chart.Correlation(sub, histogram=FALSE, method = "pearson")
dev.off()

#### Supplementary Fig. 2 SARS-CoV-2 infection status by baseline antibody levels for different coronaviruses. ####
pdf(file = "/Users/iwester4/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Antibody profiles/Manuscript/Supplementary Figure 2 SARS-CoV-2 infection status by baseline antibody levels for different coronaviruses.pdf", width = 10.5, height = 5.5)
AllParticipants_Infection %>% dplyr::select(ref_id, SEC_TRANS, contains(c("_start"))) %>% filter(!(is.na("NL63-S1-T_start") | `NL63-ecto_start`=="NA")) %>% 
  dplyr::select(-contains(c("B_", "DBS_", "MERS", "Influenza", "SARS-ecto", "SARS-S1", "SARS-CoV-2 WT RBD", "SARS-NP"))) %>% 
  group_by(ref_id) %>%
  pivot_longer(-c(ref_id, SEC_TRANS), names_to = "Virus", values_to = "value") %>% 
  mutate(SEC_TRANS = ifelse(SEC_TRANS=="Secondary transmission", "Infected", "Not infected")) %>% 
  mutate(target = ifelse(grepl("S1", Virus), "S1", ifelse(grepl("NP|np", Virus), "NP","Ecto"))) %>% 
  mutate(Virus = str_replace(Virus, "_start", "") ,Virus = str_extract(Virus, "^[^-]+")) %>% mutate(Virus = ifelse(Virus=="SARS", "SARS2", Virus)) %>% 
  mutate(Virus = factor(Virus, levels=c("SARS2", "NL63", "229E", "HKU1", "OC43"))) %>%
  ggplot(aes(x=SEC_TRANS, y=value)) + 
  scale_y_continuous(expand = c(0.1, 0), limits = c(0, 85000), breaks = seq(0, 60000, 20000)) +
  geom_jitter(aes(colour=SEC_TRANS),height = 0, width = 0.2, alpha = 0.7) + geom_boxplot(outlier.shape = NA, alpha=0.1) + theme_bw() +
  facet_grid(target ~ Virus) + scale_colour_brewer(palette = "Set1") +
  stat_compare_means(comparisons = list(c("Infected", "Not infected")), hide.ns = TRUE, label = "p.signif", label.y = 72000, tip.length = 0.01) + #, geom = "pointrange")  +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), strip.text = element_text(face = "bold")) +
  labs(y="Antibody titer", colour = "Infection status")
dev.off()

#### Supplementary Fig. 3 Disease severity by baseline antibody levels for different coronaviruses, among SARS-CoV-2 infected household members. ####
pdf(file = "/Users/iwester4/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Antibody profiles/Manuscript/Supplementary Figure 3 Disease severity by baseline antibody levels for different coronaviruses, among SARS-CoV-2 infected household members.pdf", width = 10.5, height = 5.5)
AllParticipants %>% dplyr::select(ref_id, contains("_start"), Diseaseseverity) %>% 
  dplyr::select(-contains(c("B_", "MERS", "Influenza", "SARS-ecto_start", "SARS-S1_start", "SARS-CoV-2 WT RBD_start", "SARS-NP_start"))) %>% 
  group_by(ref_id, Diseaseseverity) %>%
  pivot_longer(contains("_start"), names_to = "Virus", values_to = "value") %>% 
  mutate(target = ifelse(grepl("S1", Virus), "S1", ifelse(grepl("NP|np", Virus), "NP","Ecto"))) %>% 
  mutate(Virus = str_replace(Virus, "_start", "") ,Virus = str_extract(Virus, "^[^-]+")) %>% mutate(Virus = ifelse(Virus=="SARS", "SARS2", Virus)) %>% 
  mutate(Virus = factor(Virus, levels=c("SARS2", "NL63", "229E", "HKU1", "OC43"))) %>%
  ggplot(aes(x=Diseaseseverity, y=value)) + 
  scale_y_continuous(expand = c(0.1, 0), limits = c(0, 85000), breaks = seq(0, 60000, 20000)) +
  geom_jitter(aes(colour=Diseaseseverity),height = 0, width = 0.2, alpha = 0.7) + geom_boxplot(outlier.shape = NA, alpha=0.1) + theme_bw() +
  facet_grid(target ~ Virus)+ scale_colour_brewer(palette = "Set1", labels=c('Symptomatic', 'Pauci-symptomatic', 'Asymptomatic')) +
  stat_compare_means(comparisons = list(c("ARI", "Mild symptoms"), c("Mild symptoms", "Asymptomatic"), c("ARI", "Asymptomatic")),  hide.ns = TRUE, 
                     label = "p.signif",label.y = c(68000, 74000, 79000), geom = "pointrange")  +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), strip.text = element_text(face = "bold")) +
  labs(y="Antibody titer")
dev.off()


