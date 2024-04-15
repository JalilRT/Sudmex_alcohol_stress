#  DBM - invivo

#rm(list=ls())
# packages´ ---------------------------------------------------------------

library(tidyverse)
library(RMINC)
library(dplyr)
library(parallel)
library(lmerTest)
library(lme4)

# directories -------------------------------------------------------------

getwd()
dir.create("data", recursive = T)

# data --------------------------------------------------------------------

Jdata <- read_csv("data/DBM_dataset.csv") %>% filter( IN == "yes", Group == "Alc" | Group == "Ctrl")

#define what is the type of variables%>%

Jdata <- Jdata %>% mutate(Group = factor(Group) %>% relevel(Group, ref = "Ctrl"), # Ctrl as reference
			   Age = as.numeric(Age),
			   Sex = factor(Sex),
			   Session = factor(Session),
			   RID = factor(RID),
			   Batch = factor(Batch),
			   Intake = factor(Intake) %>% relevel(Intake , ref = "Ctrl")) # Ctrl as reference
			 
# Anatomical mask of selected ROIs

mask <- "/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo/DBM/data/ROI_mask_bin.mnc"

# Modelling

Model1 <- mincLmer(Subject ~ poly(Age,2)*Group + Sex  + Batch + (1 |RID), 
                  data = Jdata, 
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)
                  
fdrMod1 <- (mincFDR(mincLmerEstimateDF(model = Model1), mask = mask)) 

###

Model2 <- mincLmer(Subject ~ poly(Age,2)*Intake + Sex  + Batch + (1 |RID),
                  data = Jdata,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod2 <- (mincFDR(mincLmerEstimateDF(model = Model2), mask = mask))

###

Model3 <- mincLmer(Subject ~ poly(Age,2)*Group*Sex  + Batch + (1 |RID),
                  data = Jdata,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod3 <- (mincFDR(mincLmerEstimateDF(model = Model3), mask = mask))

###

Model4 <- mincLmer(Subject ~ poly(Age,2)*Intake*Sex  + Batch + (1 |RID),
                  data = Jdata,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod4 <- (mincFDR(mincLmerEstimateDF(model = Model4), mask = mask))

###

Model_ses <- mincLmer(Subject ~ Session*Group + Sex  + Batch + (1 |RID), 
                      data = Jdata, 
                      mask = mask,
                      parallel = c("local", 50),
                      REML = TRUE)
                  
fdrMod_ses <- (mincFDR(mincLmerEstimateDF(model = Model_ses), mask = mask)) 

###

Model_class_ses <- mincLmer(Subject ~ Session*Intake + Sex  + Batch + (1 |RID),
                  data = Jdata,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod_class_ses <- (mincFDR(mincLmerEstimateDF(model = Model_class_ses), mask = mask))

############### Changing reference ################

#define what is the type of variables%>%

Jdata_low <- Jdata %>% mutate(Group = factor(Group) %>% relevel(Group, ref = "Ctrl"), # Ctrl as reference
                           Age = as.numeric(Age),
                           Sex = factor(Sex),
                           Session = factor(Session),
                           RID = factor(RID),
                           Batch = factor(Batch),
                           Intake = factor(Intake) %>% relevel(Intake , ref = "Low")) # Low intake as reference

# Modelling 

Model5 <- mincLmer(Subject ~ poly(Age,2)*Intake + Sex  + Batch + (1 |RID),
                  data = Jdata_low,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod5 <- (mincFDR(mincLmerEstimateDF(model = Model5), mask = mask))

###

Model6 <- mincLmer(Subject ~ poly(Age,2)*Intake*Sex  + Batch + (1 |RID),
                  data = Jdata_low,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod6 <- (mincFDR(mincLmerEstimateDF(model = Model6), mask = mask))

# Saving data

save(Jdata, Model1, fdrMod1, Model2, fdrMod2, Model3, fdrMod3, Model4, fdrMod4, file = "DBM_data1.RData")
save(Model5, fdrMod5, Model6, fdrMod6, Model_ses, fdrMod_ses, Model_class_ses, fdrMod_class_ses, file = "DBM_data2.RData")

dir.create("tmaps/")

# Exporting

####
mincWriteVolume(Model1,"tmaps/M1-tvalue-Age_GroupAlc1.mnc", 
                column = 'tvalue-poly(Age, 2)1:GroupAlc', clobber = TRUE)

mincWriteVolume(Model1,"tmaps/M1-tvalue-Age_GroupAlc2.mnc",
                column = 'tvalue-poly(Age, 2)2:GroupAlc', clobber = TRUE)

####

mincWriteVolume(Model2,"tmaps/M2-tvalue-Age_IntakeHigh1.mnc", 
                column = 'tvalue-poly(Age, 2)1:IntakeHigh', clobber = TRUE)

mincWriteVolume(Model2,"tmaps/M2-tvalue-Age_IntakeHigh2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeHigh', clobber = TRUE)

mincWriteVolume(Model2,"tmaps/M2-tvalue-Age_IntakeLow1.mnc", 
                column = 'tvalue-poly(Age, 2)1:IntakeLow', clobber = TRUE)

mincWriteVolume(Model2,"tmaps/M2-tvalue-Age_IntakeLow2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeLow', clobber = TRUE)

####

mincWriteVolume(Model3,"tmaps/M3-tvalue-Age_GroupAlc_Sexmale1.mnc",
                column = 'tvalue-poly(Age, 2)1:GroupAlc:Sexmale', clobber = TRUE)

mincWriteVolume(Model3,"tmaps/M3-tvalue-Age_GroupAlc_Sexmale2.mnc",
                column = 'tvalue-poly(Age, 2)2:GroupAlc:Sexmale', clobber = TRUE)

####

mincWriteVolume(Model4,"tmaps/M4-tvalue-Age_IntakeHigh_Sexmale1.mnc",
                column = 'tvalue-poly(Age, 2)1:IntakeHigh:Sexmale', clobber = TRUE)

mincWriteVolume(Model4,"tmaps/M4-tvalue-Age_IntakeHigh_Sexmale2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeHigh:Sexmale', clobber = TRUE)

mincWriteVolume(Model4,"tmaps/M4-tvalue-Age_IntakeLow_Sexmale1.mnc",
                column = 'tvalue-poly(Age, 2)1:IntakeLow:Sexmale', clobber = TRUE)

mincWriteVolume(Model4,"tmaps/M4-tvalue-Age_IntakeLow_Sexmale2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeLow:Sexmale', clobber = TRUE)

####

mincWriteVolume(Model5,"tmaps/M5-tvalue-Age_IntakeHigh1.mnc",
                column = 'tvalue-poly(Age, 2)1:IntakeHigh', clobber = TRUE)

mincWriteVolume(Model5,"tmaps/M5-tvalue-Age_IntakeHigh2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeHigh', clobber = TRUE)

####

mincWriteVolume(Model6,"tmaps/M6-tvalue-Age_IntakeHigh_Sexmale1.mnc",
                column = 'tvalue-poly(Age, 2)1:IntakeHigh:Sexmale', clobber = TRUE)

mincWriteVolume(Model6,"tmaps/M6-tvalue-Age_IntakeHigh_Sexmale2.mnc",
                column = 'tvalue-poly(Age, 2)2:IntakeHigh:Sexmale', clobber = TRUE)

####

mincWriteVolume(Model_ses,"tmaps/tvalue-Sessionses-T2_GroupAlc.mnc",
                column = 'tvalue-Sessionses-T2:GroupAlc')

mincWriteVolume(Model_ses,"tmaps/tvalue-Sessionses-T3_GroupAlc.mnc",
                column = 'tvalue-Sessionses-T3:GroupAlc')

mincWriteVolume(Model_ses,"tmaps/tvalue-Sessionses-T5_GroupAlc.mnc",
                column = 'tvalue-Sessionses-T5:GroupAlc')

####
print("finished models")
