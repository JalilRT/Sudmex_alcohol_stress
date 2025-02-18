#  DBM - invivo

#rm(list=ls())
# packages´ ---------------------------------------------------------------

library(tidyverse)
library(RMINC)
library(dplyr)
library(parallel)
library(lmerTest)
library(lme4)

args <- commandArgs()
smooth <- args[6]

# directories -------------------------------------------------------------

path <- "/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo"
# Anatomical mask of selected ROIs
#mask=paste0(path,"/DBM/tomodel/atlas_labels_registered.mnc")
mask=paste0(path,"/DBM/tomodel/atlas_alllabels_registered.mnc")
setwd(paste0(path,"/DBM"))
dir.create("data", recursive = T)

# data --------------------------------------------------------------------

Jdata <- read_csv(paste0(path,"/DBM/data/DBM_dataset_",smooth,".csv")) %>% filter( IN == "yes")

#define what is the type of variables%>%

Jdata <- Jdata %>% mutate(Group = factor(Group) %>% relevel(Group, ref = "Ctrl"), # Ctrl as reference
			   Age = as.numeric(Age),
			   Sex = factor(Sex),
			   Session = factor(Session),
			   RID = factor(RID),
			   Batch = factor(Batch)) 
			 
# Modelling

Mod1 <- mincLmer(Subject ~ poly(Age,2)*Group + Sex  + Batch + (1 |RID), 
                  data = Jdata, 
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)
                  
fdrMod1 <- (mincFDR(mincLmerEstimateDF(model = Mod1), mask = mask)) 

###

Mod3 <- mincLmer(Subject ~ poly(Age,2)*Group*Sex  + Batch + (1 |RID),
                  data = Jdata,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod3 <- (mincFDR(mincLmerEstimateDF(model = Mod3), mask = mask))

###

Mod_ses <- mincLmer(Subject ~ Session*Group + Sex  + Batch + (1 |RID), 
                      data = Jdata, 
                      mask = mask,
                      parallel = c("local", 50),
                      REML = TRUE)
                  
fdrMod_ses <- (mincFDR(mincLmerEstimateDF(model = Mod_ses), mask = mask)) 

############### Changing reference ################

#define what is the type of variables%>%

Jdata_stress <- Jdata %>% mutate(Group = factor(Group) %>% relevel(Group, ref = "Str"), # Ctrl as reference
                           Age = as.numeric(Age),
                           Sex = factor(Sex),
                           Session = factor(Session),
                           RID = factor(RID),
                           Batch = factor(Batch) )

# Modelling 

Mod4 <- mincLmer(Subject ~ poly(Age,2)*Group + Sex + Batch + (1 |RID),
                  data = Jdata_stress,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod4 <- (mincFDR(mincLmerEstimateDF(model = Mod4), mask = mask))

###

Mod5 <- mincLmer(Subject ~ poly(Age,2)*Group*Sex + Batch + (1 |RID),
                  data = Jdata_stress,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod5 <- (mincFDR(mincLmerEstimateDF(model = Mod5), mask = mask))

###

Jdata_alcohol <- Jdata %>% mutate(Group = factor(Group) %>% relevel(Group, ref = "Alc"), # Ctrl as reference
                           Age = as.numeric(Age),
                           Sex = factor(Sex),
                           Session = factor(Session),
                           RID = factor(RID),
                           Batch = factor(Batch) )

Mod6 <- mincLmer(Subject ~ poly(Age,2)*Group + Sex + Batch + (1 |RID),
                  data = Jdata_alcohol,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod6 <- (mincFDR(mincLmerEstimateDF(model = Mod6), mask = mask))

###

Mod7 <- mincLmer(Subject ~ poly(Age,2)*Group*Sex + Batch + (1 |RID),
                  data = Jdata_alcohol,
                  mask = mask,
                  parallel = c("local", 50),
                  REML = TRUE)

fdrMod7 <- (mincFDR(mincLmerEstimateDF(model = Mod7), mask = mask))

# Saving data

save(Jdata, Mod1, fdrMod1, Mod3, fdrMod3, Mod4, fdrMod4, Mod5, fdrMod5, file =  paste0(path,"/DBM/smooth_",smooth,"/DBM_data1.RData"))
save(Mod6, fdrMod6, Mod7, fdrMod7, Mod_ses, fdrMod_ses, file =  paste0(path,"/DBM/smooth_",smooth,"/DBM_data2.RData"))

dir.create(paste0(path,"/DBM/smooth_",smooth,"/tmaps/"))

# Exporting

####
mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupAlc1.mnc"), 
                column = 'tvalue-poly(Age, 2)1:GroupAlc', clobber = TRUE)

mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupAlc2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc', clobber = TRUE)

mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupAlc+Str1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str', clobber = TRUE)

mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupAlc+Str2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str', clobber = TRUE)

mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupStr1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupStr', clobber = TRUE)

mincWriteVolume(Mod1,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M1-tvalue-Age_GroupStr2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupStr', clobber = TRUE)

####

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupAlc_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc:Sexmale', clobber = TRUE)

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupAlc_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc:Sexmale', clobber = TRUE)

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupAlc+Str_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str:Sexmale', clobber = TRUE)

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupAlc+Str_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str:Sexmale', clobber = TRUE)

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupStr_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupStr:Sexmale', clobber = TRUE)

mincWriteVolume(Mod3,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M3-tvalue-Age_GroupStr_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupStr:Sexmale', clobber = TRUE)

####

mincWriteVolume(Mod4,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M4-tvalue-Age_GroupAlc1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc', clobber = TRUE)

mincWriteVolume(Mod4,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M4-tvalue-Age_GroupAlc2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc', clobber = TRUE)

mincWriteVolume(Mod4,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M4-tvalue-Age_GroupAlc+Str1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str', clobber = TRUE)

mincWriteVolume(Mod4,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M4-tvalue-Age_GroupAlc+Str2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str', clobber = TRUE)

####

mincWriteVolume(Mod5,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M5-tvalue-Age_GroupAlc_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc:Sexmale', clobber = TRUE)

mincWriteVolume(Mod5,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M5-tvalue-Age_GroupAlc_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc:Sexmale', clobber = TRUE)

mincWriteVolume(Mod5,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M5-tvalue-Age_GroupAlc+Str_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str:Sexmale', clobber = TRUE)

mincWriteVolume(Mod5,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M5-tvalue-Age_GroupAlc+Str_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str:Sexmale', clobber = TRUE)

####

mincWriteVolume(Mod6,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M6-tvalue-Age_GroupAlc+Str1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str', clobber = TRUE)

mincWriteVolume(Mod6,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M6-tvalue-Age_GroupAlc+Str2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str', clobber = TRUE)

mincWriteVolume(Mod6,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M6-tvalue-Age_GroupStr1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupStr', clobber = TRUE)

mincWriteVolume(Mod6,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M6-tvalue-Age_GroupStr2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupStr', clobber = TRUE)

####

mincWriteVolume(Mod7,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M7-tvalue-Age_GroupAlc+Str_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupAlc+Str:Sexmale', clobber = TRUE)

mincWriteVolume(Mod7,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M7-tvalue-Age_GroupAlc+Str_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupAlc+Str:Sexmale', clobber = TRUE)

mincWriteVolume(Mod7,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M7-tvalue-Age_GroupStr_Sexmale1.mnc"),
                column = 'tvalue-poly(Age, 2)1:GroupStr:Sexmale', clobber = TRUE)

mincWriteVolume(Mod7,paste0(path,"/DBM/smooth_",smooth,"/tmaps/M7-tvalue-Age_GroupStr_Sexmale2.mnc"),
                column = 'tvalue-poly(Age, 2)2:GroupStr:Sexmale', clobber = TRUE)

####

#mincWriteVolume(Mod_ses,paste0(path,"/DBM/smooth_",smooth,"/tmaps/tvalue-Sessionses-T2_GroupAlc.mnc",
#                column = 'tvalue-Sessionses-T2:GroupAlc')

#mincWriteVolume(Mod_ses,paste0(path,"/DBM/smooth_",smooth,"/tmaps/tvalue-Sessionses-T3_GroupAlc.mnc",
#                column = 'tvalue-Sessionses-T3:GroupAlc')

#mincWriteVolume(Mod_ses,paste0(path,"/DBM/smooth_",smooth,"/tmaps/tvalue-Sessionses-T5_GroupAlc.mnc",
#                column = 'tvalue-Sessionses-T5:GroupAlc')

####
print("finished models")
