#  DBM - invivo

#rm(list=ls())
# packages´ ---------------------------------------------------------------

library(tidyverse)
library(RMINC)
library(dplyr)
library(parallel)
library(lmerTest)
library(lme4)
library(magrittr)
library(janitor)
library(ggpubr)
library(effects)
library(grid)
library(MRIcrotome)

# directories -------------------------------------------------------------

# Load Tmaps --------------------------------------------------------------

path="/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/"
setwd(path)
anatVol <- mincArray(mincGetVolume(paste0(path,"Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template_Masked.mnc")))
#anatVol <- mincArray(mincGetVolume(paste0(path,"analysis/smri/DBM_invivo/DBM/data/secondlevel_template0_maskBrainExtractionBrain.mnc")))

# Plot Tmaps --------------------------------------------------------------

dir.create(paste0(path,"analysis/fmri/rabies/melodic/fig_tmaps"))

############## 1 slice only significant models ############

nets=c("0009","0010")

for (Net in nets) {

  melodic <- mincArray(mincGetVolume(paste0(path,"analysis/fmri/rabies/melodic/split/melodic_",Net,".mnc")) )

  svg(paste0("analysis/fmri/rabies/melodic/fig_tmaps/single/melodic_",Net,"_sag.svg"), height = 3, width = 3, bg = "transparent")
    sliceSeries(nrow = 1, ncol = 1, begin = 62, end = 62, dimension = 1) %>%
      anatomy(anatVol, low=2000, high=20000) %>%
      overlay(melodic, low=2, high=33, symmetric = F) %>% 
    draw()
  dev.off()

  svg(paste0("analysis/fmri/rabies/melodic/fig_tmaps/single/melodic_",Net,"_cor.svg"), height = 3, width = 3, bg = "transparent")
    sliceSeries(nrow = 1, ncol = 1, begin = 80, end = 80) %>%
      anatomy(anatVol, low=2000, high=20000) %>%
      overlay(melodic, low=2, high=33, symmetric = T) %>% 
    draw()
  dev.off()

  svg(paste0("analysis/fmri/rabies/melodic/fig_tmaps/single/melodic_",Net,"_cor2.svg"), height = 3, width = 3, bg = "transparent")
    sliceSeries(nrow = 1, ncol = 1, begin = 130, end = 130) %>%
      anatomy(anatVol, low=2000, high=20000) %>%
      overlay(melodic, low=2, high=33, symmetric = T) %>% 
    draw()
  dev.off()

  svg(paste0("analysis/fmri/rabies/melodic/fig_tmaps/single/melodic_",Net,"_axial.svg"), height = 3, width = 3, bg = "transparent")
    sliceSeries(nrow = 1, ncol = 1, begin = 80, end = 80, dimension = 3) %>%
      anatomy(anatVol, low=0.5, high=4) %>%
      overlay(melodic, low=2, high=33, symmetric = T) %>% 
    draw()
  dev.off()


}
