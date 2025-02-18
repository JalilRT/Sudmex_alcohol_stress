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

args <- commandArgs()
smooth <- args[6]

# directories -------------------------------------------------------------
smooth="1mm"
path <- "/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo"
atlas=paste0(path,"/DBM/tomodel/atlas_alllabels_registered.mnc")
mask=paste0(path,"/jacobians/output/secondlevel/final/average/mask_shapeupdate.mnc")
anatVol=mincArray(mincGetVolume(paste0(path,"/DBM/tomodel/template_sharpen_shapeupdate_brain.mnc")))
setwd(paste0(path,"/DBM/mask_less/smooth_",smooth))
dir.create("data", recursive = T)

# Load data ---------------------------------------------------------------

load("DBM_data1.RData")
load("DBM_data2.RData")

# Some specifications -----------------------------------------------------

dir.create("fig_tmaps")
system("rm fig_tmaps/*tsv")
system("for i in $(ls $PWD/tmaps/M*mnc); do echo $(basename $i .mnc); done > fig_tmaps/names_slices.tsv")
system("for i in $(ls $PWD/tmaps/M*mnc); do mincbbox -minccrop $i; done > fig_tmaps/slices.tsv")

names_slices <- read_table("fig_tmaps/names_slices.tsv", col_names=F)
slices <- read_table("fig_tmaps/slices.tsv", col_names=F) %>% add_column(names_slices %>% rename("N1" = "X1"), .before=1)

# Plot Tmaps --------------------------------------------------------------

Mods=c("Mod1","Mod3","Mod4","Mod5","Mod6","Mod7") #modify accordig the models
Groups=c("Alc","Str","Alc+Str") #modify accordig the groups
pal_group=list(c("#C96505"),c("#00A087"),c("#09467C")) %>% set_names(Groups)
  #set_names(c("Ctrl","Alc","Str","Alc+Str"))c("#757575"),

Allthresholds <- Mods %>% map(~ paste0("fdr",.x) ) %>% map(~ attr(get(.x), "thresholds") ) %>% set_names(Mods)
AlldimThr <- Allthresholds %>% map(~ dimnames(.x)[[2]][] )# %>% set_names(Mods)

Mods_wS <- Mods[which(map_lgl(AlldimThr, ~ {any(str_detect(.x,":Sex")) }) == "FALSE")]
Mods_S <- Mods[which(map_lgl(AlldimThr, ~ {any(str_detect(.x,":Sex")) }) == "TRUE")]

# Plot Tmaps --------------------------------------------------------------

# For cycle
for (Mod in Mods_wS) {

 rois=read_csv(paste0(path,"/DBM/smooth_",smooth,"/Trayectories/ROIs.csv")) %>%
  filter(grepl(paste0("L",Mod),Model))

 fdrMod=paste0("fdr",Mod)
 thresholds=attr(get(fdrMod), "thresholds")
 dimThr=dimnames(thresholds)[[2]][c(-1)]

for (Group in Groups) {

if (paste0("tvalue-poly(Age, 2)1:Group",Group) %in% dimThr) {

  if (Group == Groups[1]) {
    
    if (paste0("tvalue-poly(Age, 2)1:Group",Groups[1]) %in% dimThr) {
    assign(paste0(Mod, "_stats1_", Groups[1]), mincArray(get(Mod), paste0("tvalue-poly(Age, 2)1:Group",Groups[1])))
    } else {
      print("No contrast for this group")
    }

  } else if (Group == Groups[2]) {

    if (paste0("tvalue-poly(Age, 2)1:Group",Groups[2]) %in% dimThr) {
    assign(paste0(Mod, "_stats1_", Groups[2]), mincArray(get(Mod), paste0("tvalue-poly(Age, 2)1:Group",Groups[2])))
    } else {
      print("No contrast for this group")
    }

  } else if (Group == Groups[3]) {

    if (paste0("tvalue-poly(Age, 2)1:Group",Groups[3]) %in% dimThr) {
    assign(paste0(Mod, "_stats1_", Groups[3]), mincArray(get(Mod), paste0("tvalue-poly(Age, 2)1:Group",Groups[3])))
    } else {
      print("No contrast for this group")
    }

  }

  # ## Checking groups and stats

  # if (Group == Groups[1]) {
    
  #   if (paste0(Mod,"_stats1_",Groups[1]) %in% ls()) {
  #     stats1 <- get(paste0(Mod,"_stats1_",Groups[1]))
  #   } else {
  #     print("No contrast for this group")
  #   }

  # } else if (Group == Groups[2]) {

  #   if (paste0(Mod,"_stats1_",Groups[2]) %in% ls()) {
  #     stats1 <- get(paste0(Mod,"_stats1_",Groups[2]))
  #   } else {
  #     print("No contrast for this group")
  #   }

  # } else if (Group == Groups[3]) {

  #   if (paste0(Mod,"_stats1_",Groups[3]) %in% ls()) {
  #     stats1 <- get(paste0(Mod,"_stats1_",Groups[3]))
  #   } else {
  #     print("No contrast for this group")
  #   }
  # }
  
  tFDR01L <- Allthresholds[[Mod]]["0.01", paste0("tvalue-poly(Age, 2)1:Group",Group)] %>% round(digits=2)
  tFDR5L <- Allthresholds[[Mod]]["0.05", paste0("tvalue-poly(Age, 2)1:Group",Group)] %>% round(digits=2)
  tFDR1L <- Allthresholds[[Mod]]["0.1", paste0("tvalue-poly(Age, 2)1:Group",Group)] %>% round(digits=2)
  tFDR2L <- Allthresholds[[Mod]]["0.2", paste0("tvalue-poly(Age, 2)1:Group",Group)] %>% round(digits=2)
  tmax <- 7
  
  dir.create(paste0("fig_tmaps/single/",Mod,"/"), recursive = T)

  if (tFDR2L %in% NA){
    print("non significant at 20% FDR")
  } else {
    dir.create("fig_tmaps/single",recursive=T)
    for (roi in 1:nrow(rois)) {
      svg(paste0("fig_tmaps/single/",Mod,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
      sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
      anatomy(anatVol, low=0.15, high=3) %>%   
      overlay(get(paste0(Mod, "_stats1_", Group) ), 
        low=tFDR2L, high=tmax, symmetric = T) %>%
      draw()
      dev.off()
    }
  }

  if (tFDR5L %in% NA){
    print("non significant at 5% FDR")
  } else {
    dir.create("fig_tmaps/single",recursive=T)
    for (roi in 1:nrow(rois)) {
      svg(paste0("fig_tmaps/single/",Mod,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
      sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
      anatomy(anatVol, low=0.15, high=3) %>%   
      overlay(get(paste0(Mod, "_stats1_", Group) ),
        #col = colorRampPalette(c("#5CB270", "yellow"))(100),rCol = colorRampPalette(c("#6407B0", "yellow"))(100), 
        low=tFDR2L, high=tmax, symmetric = T) %>%
      contours(get(paste0(Mod, "_stats1_", Group) ), 
        levels=tFDR5L, col="yellow", lty=5) %>% 
      draw()
      dev.off()
    }
  }

  if (tFDR01L %in% NA){
    print("non significant at 1% FDR")
  } else {
    dir.create("fig_tmaps/single",recursive=T)
    for (roi in 1:nrow(rois)) {
      svg(paste0("fig_tmaps/single/",Mod,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
      sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
      anatomy(anatVol, low=0.15, high=3) %>%   
      overlay(get(paste0(Mod, "_stats1_", Group) ), low=tFDR2L, high=tmax, symmetric = T) %>%
      contours(get(paste0(Mod, "_stats1_", Group) ), levels=c(tFDR5L,tFDR01L), col=c("yellow","green"), lty=5) %>%
      draw()
      dev.off()
    }
  }

} else {
  print("No contrast for this group")
}

} }

# For cycle
for (Mod in Mods_S) {

for (Group in Groups) {

  rois=read_csv(paste0(path,"/DBM/mask_less/smooth_",smooth,"/Trayectories/ROIs.csv")) %>%
  filter(grepl(paste0("L",Mod),Model)) 

  fdrMod=paste0("fdr",Mod)
  thresholds=attr(get(fdrMod), "thresholds")
  dimThr=dimnames(thresholds)[[2]][c(-1)]

if (paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale") %in% dimThr) {
 
  if (grepl("\\+", Group)) {
    filter_group <- gsub("\\+", "", Group)
  } else {
    filter_group <- Group
  }

  rois <- rois %>%
    filter(grepl(paste0("_",filter_group,"$"), rois$Model))

  print(paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale"))  
  assign(paste0(Mod, "_stats1_", Group), mincArray(get(Mod), paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale")))

  tFDR01L <- Allthresholds[[Mod]]["0.01", paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale")] %>% round(digits=2)
  tFDR5L <- Allthresholds[[Mod]]["0.05", paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale")] %>% round(digits=2)
  tFDR1L <- Allthresholds[[Mod]]["0.1", paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale")] %>% round(digits=2)
  tFDR2L <- Allthresholds[[Mod]]["0.2", paste0("tvalue-poly(Age, 2)1:Group",Group,":Sexmale")] %>% round(digits=2)
  tmax <- 7
  
  

  # if (tFDR2L %in% NA){
  #   print("non significant at 20% FDR")
  # } else {
  #   dir.create(paste0("fig_tmaps/single/",Mod,"/",Group), recursive = T)
  #   for (roi in 1:nrow(rois)) {
  #     svg(paste0("fig_tmaps/single/",Mod,"/",Group,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
  #     sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
  #     anatomy(anatVol, low=0.15, high=3) %>%   
  #     overlay(get(paste0(Mod, "_stats1_", Group)),
  #       low=tFDR2L, high=tmax, symmetric = T) %>%
  #     draw()
  #     dev.off()
  #   }
  # }

  if (tFDR5L %in% NA){
    print("non significant at 5% FDR")
  } else {
    dir.create(paste0("fig_tmaps/single/",Mod,"/",Group), recursive = T)
    for (roi in 1:nrow(rois)) {
      svg(paste0("fig_tmaps/single/",Mod,"/",Group,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
      sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
      anatomy(anatVol, low=0.15, high=3) %>%   
      overlay(get(paste0(Mod, "_stats1_", Group)), 
        low=tFDR2L, high=tmax, symmetric = T) %>%
      contours(abs(get(paste0(Mod, "_stats1_", Group))), 
        levels=tFDR5L, col="yellow", lty=6) %>% 
      draw()
      dev.off()
    }
  }

  if (tFDR01L %in% NA){
    print("non significant at 1% FDR")
  } else {
    dir.create(paste0("fig_tmaps/single/",Mod,"/",Group), recursive = T)
    for (roi in 1:nrow(rois)) {
      svg(paste0("fig_tmaps/single/",Mod,"/",Group,"/",Mod,"-",rois$ROI_model[roi],".svg"), height = 3, width = 3, bg = "transparent")
      sliceSeries(nrow = 1, ncol = 1, begin = rois$d2[roi], end = rois$d2[roi]) %>%
      anatomy(anatVol, low=0.15, high=3) %>%   
      overlay(get(paste0(Mod, "_stats1_", Group)), 
        low=tFDR2L, high=tmax, symmetric = T) %>%
      contours(abs(get(paste0(Mod, "_stats1_", Group))), 
        levels=c(tFDR5L,tFDR01L), col=c("yellow","green"), lty=6) %>%
      draw()
      dev.off()
    }
  }

} else {
  print("No contrast for this group")
}

} }


