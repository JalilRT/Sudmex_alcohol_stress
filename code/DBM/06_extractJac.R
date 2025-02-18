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

args <- commandArgs()
smooth <- args[6]

# directories -------------------------------------------------------------
#rm(list = ls())
path <- "/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo"
atlas=paste0(path,"/DBM/tomodel/atlas_alllabels_registered.mnc")
mask=paste0(path,"/jacobians/output/secondlevel/final/average/mask_shapeupdate.mnc")
setwd(paste0(path,"/DBM/smooth_",smooth))
dir.create("data", recursive = T)

# Load data ---------------------------------------------------------------

load("DBM_data1.RData")
load("DBM_data2.RData")

# Extract peaks
Mods=c("Mod1","Mod3","Mod5","Mod7","Mod4","Mod6") #modify accordig the models
Groups=c("Alc","Str","Alc+Str") #
Allthresholds <- Mods %>% map(~ paste0("fdr",.x) ) %>% map(~ attr(get(.x), "thresholds") ) %>% set_names(Mods)
dimThr <- Allthresholds %>% map(~ dimnames(.x)[[2]][] )# %>% set_names(Mods)

Mods_wS <- Mods[which(map_lgl(dimThr, ~ {any(str_detect(.x,":Sex")) }) == "FALSE")]
Mods_S <- Mods[which(map_lgl(dimThr, ~ {any(str_detect(.x,":Sex")) }) == "TRUE")]

## For models without Sex interaction (Age*Group + covs)
predictorL_wS <- Groups %>% map(function(y) dimThr %>% map(~ match(paste0("tvalue-poly(Age, 2)1:Group",y),.x) )) %>% set_names(Groups)
columnL_wS <- Groups %>% map(function(y) Mods_wS %>% map(~ dimThr[[.x]][predictorL_wS[[y]][[.x]]] ) %>% set_names(Mods_wS) ) %>% set_names(Groups) %>% map(~ .x %>% discard(is.na) ) %>% compact()
tFDR5L_wS <- Groups %>% map(function(y) Mods_wS %>% map(~ Allthresholds[[.x]]["0.05",predictorL_wS[[y]][[.x]]] ) %>% set_names(Mods_wS) ) %>% set_names(Groups) %>% map(~ .x %>% discard(is.na) ) %>% compact()
temp_wS <- names(tFDR5L_wS) %>% map(function(y) names(tFDR5L_wS[[y]]) %>% map(~ mincFindPeaks(get(.x), column = columnL_wS[[y]][[.x]], direction = "both", threshold = tFDR5L_wS[[y]][[.x]], minDistance = 4)) %>% 
	    set_names(paste0("L",names(tFDR5L_wS[[y]]))) ) %>% set_names(names(tFDR5L_wS)) %>% unlist(recursive=FALSE) 
mod_list_wS <- temp_wS %>% set_names(paste0(names(temp_wS) %>% str_split(pattern='\\.', n=2) %>% 
	    map_chr(., 2),"_",names(temp_wS) %>% str_split(pattern='\\.', n=2) %>% map_chr(., 1)) %>% str_replace("\\+","") )

# ## For models with Sex interaction (Age*Group*Sex + covs)
predictorL_S <- Groups %>% map(function(y) dimThr %>% map(~ match(paste0("tvalue-poly(Age, 2)1:Group",y,":Sexmale"),.x) )) %>% set_names(Groups)
columnL_S <- Groups %>% map(function(y) Mods_S %>% map(~ dimThr[[.x]][predictorL_S[[y]][[.x]]] ) %>% set_names(Mods_S) ) %>% set_names(Groups)
tFDR5L_S <- Groups %>% map(function(y) Mods_S %>% map(~ Allthresholds[[.x]]["0.05",predictorL_S[[y]][[.x]]] ) %>% set_names(Mods_S) ) %>% set_names(Groups) %>% map(~ .x %>% discard(is.na) ) %>% compact()
temp_S <- names(tFDR5L_S) %>% map(function(y) names(tFDR5L_S[[y]]) %>% map(~ mincFindPeaks(get(.x), column = columnL_S[[y]][[.x]], direction = "both", threshold = tFDR5L_S[[y]][[.x]], minDistance = 4)) %>% 
 	    set_names(paste0("L",names(tFDR5L_S[[y]]))) ) %>% set_names(names(tFDR5L_S)) %>% unlist(recursive=FALSE) 
mod_list_S <- temp_S %>% set_names(paste0(names(temp_S) %>% str_split(pattern='\\.', n=2) %>% 
	    map_chr(., 2),"_",names(temp_S) %>% str_split(pattern='\\.', n=2) %>% map_chr(., 1)) %>% str_replace("\\+","") )

mod_peaks <- c(mod_list_wS,mod_list_S) %>%
  map(~ mincLabelPeaks(.x,  atlas, defs="../data/SIGMA_InVivo_Anatomical_Brain_Atlas_ListOfStructures.csv") )


dir.create("Peaks")

mod_peaks %>% iwalk(~write_csv(.x, paste0(getwd(),"/Peaks/",.y, ".csv")))

# Trayectories -------------------------------------------------------------

mod_peaks <- mod_peaks %>% map(~ .x %>% mutate(label_clean = label %>% make_clean_names()) )

ROIs <- names(mod_peaks) %>% map(~ mod_peaks[[.x]] %>% add_column(Model = rep(.x, nrow(mod_peaks[[.x]]) )) ) %>% reduce(rbind) %>% mutate(ROI_model = str_c(label_clean,'_',Model) )

Jdata_jacobians <- Jdata %>% select(RID,Session,Subject,Group,Age,Sex,Batch)

for (i in seq(1:length(ROIs$x)) ) {
Jdata_jacobians[ROIs$ROI_model[i]] = mincGetWorldVoxel(Jdata$Subject, ROIs$x[i], ROIs$y[i], ROIs$z[i])
}

dir.create("Trayectories")

write_csv(Jdata_jacobians,"Trayectories/Jdata_jacobians.csv")
write_csv(ROIs,"Peaks/ROIs.csv")


# Subgrouping -------------------------------------------------------------

# By sex

Jdata_jacobians <- Jdata_jacobians %>% mutate(Group_sex = case_when(Group == "Alc" & Sex == "female" ~ "Alc F",
						Group == "Alc" & Sex == "male" ~ "Alc M",
						Group == "Ctrl" & Sex == "female" ~ "Ctrl F",
                                                Group == "Ctrl" & Sex == "male" ~ "Ctrl M",
						Group == "Str" & Sex == "female" ~ "Str F",
						Group == "Str" & Sex == "male" ~ "Str M",
						Group == "Alc+Str" & Sex == "female" ~ "Alc+Str F",
                                                Group == "Alc+Str" & Sex == "male" ~ "Alc+Str M"), .before = 8) %>%
			  mutate(Group_sex = factor(Group_sex, levels = c("Ctrl F","Ctrl M","Alc F","Alc M","Str F","Str M","Alc+Str F","Alc+Str M")) )

# Plotting ----------------------------------------------------------------

theme_settings <- theme(text = element_text(size=20),
                        axis.text.x = element_text(size=15))

# colors
#pal_group <- c(alpha("#000000",0.6),"#B22222","#217175")
pal_sex <- c(alpha("#83458E",1),alpha("#28799f",1))
pal_class <- c("#DF9A57","#B22222","#217175","#44355B")
pal_group <- c("#737373","#B22222","#168a44","#217175")
pal_groupSex <- c("#737373","#424242","#7DA2E8","#00688B","#9780F4","#6C3F81","#B22222","#217175")
                        
nROI <- ROIs$ROI_model

# Linear

dir.create("Trayectories/Linear/Group", recursive = TRUE)

plots_ROIs_L <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_L[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Group",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Group", "Age"),
                                          lmer(get(ROI) ~ Age*Group + 
                                                 Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Group), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Group", "Age"),
                                            lmer(get(ROI) ~ Age*Group + 
                                                   Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                            xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                                 max(Jdata_jacobians$Age),1)))), 
		  aes(y=fit, ymin=lower, ymax=upper, fill = Group), alpha=0.1) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}
  
1:length(plots_ROIs_L) %>% map(~ ggsave(filename = paste0("Trayectories/Linear/Group/",plots_ROIs_L[.x] %>% names(),".png"), plot = plots_ROIs_L[[.x]],dpi = 300, width = 5.5, height = 4.5) ) 


dir.create("Trayectories/Linear/Group_sex", recursive = TRUE)

plots_ROIs_sex_L <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_sex_L[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Group_sex",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group_sex", linewidth =2,
                  palette = pal_groupSex,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) + 
        geom_line(data = as_tibble(Effect(c("Group_sex","Age"),
                                          lmer(get(ROI) ~ Age*Group_sex + Batch + (1 |RID), data = Jdata_jacobians),
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))),
                  aes(y=fit, color = Group_sex), size=2) +
	labs(fill = "Group by sex", color = "Group by sex") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
        theme_settings
}

1:length(plots_ROIs_sex_L) %>% map(~ ggsave(filename = paste0("Trayectories/Linear/Group_sex/",plots_ROIs_sex_L[.x] %>% names(),".png"), plot = plots_ROIs_sex_L[[.x]],dpi = 300, width = 5.5, height = 4.5) )

#-------------------------------------------------------------------------

# Polynomial

dir.create("Trayectories/Polynomial/Group", recursive = TRUE)

plots_ROIs <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Group",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Group", "Age"),
                                          lmer(get(ROI) ~ poly(Age,2)*Group + 
                                                 Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Group), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Group", "Age"),
                                            lmer(get(ROI) ~ poly(Age,2)*Group + 
                                                   Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                            xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                                 max(Jdata_jacobians$Age),1)))), 
		  aes(y=fit, ymin=lower, ymax=upper, fill = Group), alpha=0.3) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}
  
1:length(plots_ROIs) %>% map(~ ggsave(filename = paste0("Trayectories/Polynomial/Group/",plots_ROIs[.x] %>% names(),".png"), plot = plots_ROIs[[.x]],dpi = 300, width = 5.5, height = 4.5) ) 


dir.create("Trayectories/Polynomial/Group_sex", recursive = TRUE)

plots_ROIs_sex <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_sex[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Group_sex",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group_sex", linewidth =2,
                  palette = pal_groupSex,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) + 
        geom_line(data = as_tibble(Effect(c("Group_sex","Age"),
                                          lmer(get(ROI) ~ poly(Age,2)*Group_sex + Batch + (1 |RID), data = Jdata_jacobians),
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))),
                  aes(y=fit, color = Group_sex), size=2) +
	labs(fill = "Group by sex", color = "Group by sex") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
        theme_settings
}

1:length(plots_ROIs_sex) %>% map(~ ggsave(filename = paste0("Trayectories/Polynomial/Group_sex/",plots_ROIs_sex[.x] %>% names(),".png"), plot = plots_ROIs_sex[[.x]],dpi = 300, width = 5.5, height = 4.5) )


# Creating table of peaks -------------------------------------------------

ROI_table <- ROIs %>% mutate(ROIs = gsub('_left','',label %>% 
                        str_split(pattern = " ", n = 2) %>% 
                        reduce(rbind) %>% .[,2]),
                      Hemisphere = label %>% 
                        str_split(pattern = " ", n = 2) %>% 
                        reduce(rbind) %>% .[,1],
                      "t-value" = value,
                      Coordinates = str_c(round(x,2),',',round(y,2),',',round(z,2)),
                      Effect = case_when(value > 0 ~ "Increase",
                                         TRUE ~ "Decrease"),
                      Model2 = case_when(grepl("LMod1",Model) ~ "Model 1",
                                         grepl("LMod2",Model) ~ "Model 2",
                                         grepl("LMod3",Model) ~ "Model 3",
                                         grepl("LMod4",Model) ~ "Model 4",
					 grepl("LMod5",Model) ~ "Model 5",
                                         grepl("LMod6",Model) ~ "Model 6",
                                         grepl("LMod7",Model) ~ "Model 7"),
                      Contrast = case_when(Model == "LMod1_Alc" ~ "EtOH > Ctrl",
                                        Model == "LMod1_AlcStr" ~ "EtOH+Str > Ctrl",
                                        Model == "LMod1_Str" ~ "Str > Ctrl",
                                        Model == "LMod3_Alc" ~ "EtOH > Ctrl",
                                        Model == "LMod3_Str" ~ "Str > Ctrl",
                                        Model == "LMod3_AlcStr" ~ "EtOH+Str > Ctrl",
                                        Model == "LMod4_Alc" ~ "EtOH > Str",
                                        Model == "LMod5_Alc" ~ "EtOH > Str",
                                        Model == "LMod5_AlcStr" ~ "EtOH+Str > Str",
					Model == "LMod6_AlcStr" ~ "EtOH+Str > Alc",
                                        Model == "LMod6_Str" ~ "Str > Alc",
                                        Model == "LMod7_AlcStr" ~ "EtOH+Str > Alc",
                                        Model == "LMod7_Str" ~ "Str > Alc"),
 		      Type = case_when(grepl("LMod",Model) ~ "Linear",
                                       grepl("PMod",Model) ~ "Poly") ) %>% 
  select(Model2,Type,Contrast,ROIs, Hemisphere, Coordinates, Effect, "t-value")

write_csv(ROIs,"Trayectories/ROIs.csv")
write_csv(ROI_table,"Trayectories/ROI_table.csv")

save(atlas,mod_peaks,theme_settings,plots_ROIs_L,plots_ROIs_sex_L,Jdata_jacobians,ROIs,ROI_table,nROI, file = "Trayectories_data1.RData")
save(plots_ROIs,plots_ROIs_sex, file = "Trayectories_data2.RData")
