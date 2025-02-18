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
#smooth="1mm"
path <- "/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo"
atlas=paste0(path,"/DBM/tomodel/atlas_alllabels_registered.mnc")
mask=paste0(path,"/jacobians/output/secondlevel/final/average/mask_shapeupdate.mnc")
anatVol=mincArray(mincGetVolume(paste0(path,"/DBM/tomodel/template_sharpen_shapeupdate_brain.mnc")))
setwd(paste0(path,"/DBM/smooth_",smooth))
dir.create("data", recursive = T)
#atlas = "/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.mnc"
listROI = paste0(path,"/DBM/data/SIGMA_InVivo_Anatomical_Brain_Atlas_ListOfStructures.csv")

# Load data ---------------------------------------------------------------

load("DBM_data1.RData")
load("DBM_data2.RData")
#Jdata_jacobians <- read_csv(paste0(path,"/DBM/smooth_",smooth,"/Trayectories/Jdata_jacobians.csv")) 

# Extract jacobians all ROIs ----------------------------------------------

Volumes_jacobian_ob <- anatGetAll(filenames=Jdata$Subject, 
                      atlas=atlas, 
                      method="jacobians",
                      defs=listROI)
                      
Volumes_jacobian <- Volumes_jacobian_ob %>% clean_names() %>% 
    as_tibble %>% 
    mutate(Subject = Jdata$Subject,
    RID = Jdata$RID,
    Session = Jdata$Session,
    Age = Jdata$Age,
    Group = Jdata$Group,
    Sex = Jdata$Sex,
    Batch = Jdata$Batch) %>% 
    select(RID,Group,Session,Sex,Batch,Subject,everything())

write_csv(Volumes_jacobian, "Trayectories/Volumes_jacobians.csv")

# Plot all ROIs -----------------------------------------------------------


Volumes_jacobian <- Volumes_jacobian %>% mutate(Group_sex = case_when(Group == "Alc" & Sex == "female" ~ "Alc F",
						Group == "Alc" & Sex == "male" ~ "Alc M",
						Group == "Ctrl" & Sex == "female" ~ "Ctrl F",
                                                Group == "Ctrl" & Sex == "male" ~ "Ctrl M",
						Group == "Str" & Sex == "female" ~ "Str F",
						Group == "Str" & Sex == "male" ~ "Str M",
						Group == "Alc+Str" & Sex == "female" ~ "Alc+Str F",
                                                Group == "Alc+Str" & Sex == "male" ~ "Alc+Str M"), .before = 7) %>%
			  mutate(Group_sex = factor(Group_sex, levels = c("Ctrl F","Ctrl M","Alc F","Alc M","Str F","Str M","Alc+Str F","Alc+Str M")) )

# Plotting ----------------------------------------------------------------

theme_settings <- theme(text = element_text(size=20),
                        axis.text.x = element_text(size=15))

# colors
pal_sex <- c(alpha("#83458E",1),alpha("#28799f",1))
pal_class <- c("#DF9A57","#B22222","#217175","#44355B")
pal_group <- c("#737373","#B22222","#168a44","#217175")
pal_groupSex <- c("#737373","#424242","#7DA2E8","#00688B","#9780F4","#6C3F81","#B22222","#217175")
                        
ROIs <- Volumes_jacobian %>% select( - c(RID,Group,Session,Age,Sex,Batch,Subject,Group_sex)) %>% colnames()
nROI <- ROIs

dir.create("Trayectories/ROI/Group", recursive = TRUE)

plots_ROIs_V <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_V[[ROI]] <- ggscatter(Volumes_jacobian,
                  x = "Age", y = ROI, group = "Group",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Group", "Age"),
                                          lmer(get(ROI) ~ Age*Group + 
                                                 Sex + Batch + (1 |RID), data = Volumes_jacobian), 
                                          xlevels=list(Age=seq(min(Volumes_jacobian$Age),
                                                               max(Volumes_jacobian$Age),1)))), 
                  aes(y=fit, color = Group), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Group", "Age"),
                                            lmer(get(ROI) ~ Age*Group + 
                                                   Sex + Batch + (1 |RID), data = Volumes_jacobian), 
                                            xlevels=list(Age=seq(min(Volumes_jacobian$Age),
                                                                 max(Volumes_jacobian$Age),1)))), 
		  aes(y=fit, ymin=lower, ymax=upper, fill = Group), alpha=0.1) +
        theme(#legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_V) %>% map(~ ggsave(filename = paste0("Trayectories/ROI/Group/",plots_ROIs_V[.x] %>% names(),".png"), plot = plots_ROIs_V[[.x]],dpi = 300, width = 5.5, height = 4.5) ) 

dir.create("Trayectories/ROI/Group_sex", recursive = TRUE)

plots_ROIs_V_sex <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_V_sex[[ROI]] <- ggscatter(Volumes_jacobian,
                  x = "Age", y = ROI, group = "Group",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Group", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) + facet_grid(. ~ Sex) +
        geom_line(data = as_tibble(Effect(c("Group", "Age"),
                                          lmer(get(ROI) ~ Age*Group + 
                                                 Sex + Batch + (1 |RID), data = Volumes_jacobian), 
                                          xlevels=list(Age=seq(min(Volumes_jacobian$Age),
                                                               max(Volumes_jacobian$Age),1)))), 
                  aes(y=fit, color = Group), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Group", "Age"),
                                            lmer(get(ROI) ~ Age*Group + 
                                                   Sex + Batch + (1 |RID), data = Volumes_jacobian), 
                                            xlevels=list(Age=seq(min(Volumes_jacobian$Age),
                                                                 max(Volumes_jacobian$Age),1)))), 
		  aes(y=fit, ymin=lower, ymax=upper, fill = Group), alpha=0.1) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_V_sex) %>% map(~ ggsave(filename = paste0("Trayectories/ROI/Group_sex/",plots_ROIs_V_sex[.x] %>% names(),".png"), plot = plots_ROIs_V_sex[[.x]],dpi = 300, width = 5.5, height = 4.5) ) 


# Boxplot ----------------------------------------------------------------
# my_comparisons <- list(c("Ctrl","Alc"),c("Ctrl","Str"),c("Ctrl","Alc+Str"), c("Alc","Str"),c("Alc","Alc+Str"),c("Str","Alc+Str"))

# plots_ROIs_V <- NULL
# for (i in 1:length(nROI)) {
#   ROI <- nROI[i]
#   plots_ROIs_V[[ROI]] <- ggboxplot(Volumes_jacobian, 
#   		  x = "Group", y = ROI,
# 	          color = "Group", linewidth =2,
#                   palette = pal_group,
#                   add.params = list(size = 2, alpha = 0.5),
#                   plot_type ="b", add = "jitter",
#                   title = gsub('_left','',ROIs[i]),
#                   xlab = "Group",
#                   ylab = "Local volume",
#                   font.x = c(16,"bold"),
#                   font.y = c(16,"bold"),
#                   font.tickslab = c(14,"bold"),
#                   ref.group = "Ctrl") + facet_grid(. ~ Session) +
#               theme(legend.position = "none",axis.title.x=element_blank(),
#               plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
#               theme_settings +
#   stat_compare_means(comparisons = my_comparisons, ref.group = "Ctrl") 
# }


# dir.create("Trayectories/boxplot/Group_sex/male", recursive = TRUE)

# plots_ROIs_male_V <- NULL
# for (i in 1:length(nROI)) {
#   ROI <- nROI[i]
#   plots_ROIs_sex_V[[ROI]] <- ggboxplot(Volumes_jacobian, 
#   		  x = "Group", y = ROI,
# 	          color = "Group", linewidth =2,
#                   palette = pal_group,
#                   add.params = list(size = 2, alpha = 0.5),
#                   plot_type ="b", add = "jitter",
#                   title = gsub('_left','',ROIs[i]),
#                   xlab = "Group",
#                   ylab = "Local volume",
#                   font.x = c(16,"bold"),
#                   font.y = c(16,"bold"),
#                   font.tickslab = c(14,"bold"),
#                   ref.group = "Ctrl") + facet_grid(. ~ Session) +
#               theme(legend.position = "none",axis.title.x=element_blank(),
#               plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
#               theme_settings +
#   stat_compare_means(comparisons = my_comparisons, ref.group = "Ctrl") 
# }

# 1:length(plots_ROIs_male_V) %>% map(~ ggsave(filename = paste0("Trayectories/boxplot/Group_sex/male",plots_ROIs_male_V[.x] %>% names(),".png"), plot = plots_ROIs_male_V[[.x]],dpi = 300, width = 5.5, height = 4.5) )

# plots_ROIs_female_V <- NULL
# for (i in 1:length(nROI)) {
#   ROI <- nROI[i]
#   plots_ROIs_sex_V[[ROI]] <- ggboxplot(Volumes_jacobian, 
#   		  x = "Group", y = ROI,
# 	          color = "Group", linewidth =2,
#                   palette = pal_group,
#                   add.params = list(size = 2, alpha = 0.5),
#                   plot_type ="b", add = "jitter",
#                   title = gsub('_left','',ROIs[i]),
#                   xlab = "Group",
#                   ylab = "Local volume",
#                   font.x = c(16,"bold"),
#                   font.y = c(16,"bold"),
#                   font.tickslab = c(14,"bold"),
#                   ref.group = "Ctrl") + facet_grid(. ~ Session) +
#               theme(legend.position = "none",axis.title.x=element_blank(),
#               plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) +
#               theme_settings +
#   stat_compare_means(comparisons = my_comparisons, ref.group = "Ctrl") 
# }

# 1:length(plots_ROIs_female_V) %>% map(~ ggsave(filename = paste0("Trayectories/boxplot/Group_sex/female",plots_ROIs_female_V[.x] %>% names(),".png"), plot = plots_ROIs_female_V[[.x]],dpi = 300, width = 5.5, height = 4.5) )
