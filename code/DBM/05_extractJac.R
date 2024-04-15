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

# directories -------------------------------------------------------------

getwd()

# Load data ---------------------------------------------------------------

load("DBM_data1.RData")
load("DBM_data2.RData")

# Extract peaks

atlas='/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo/DBM/data/ROI_mask.mnc'

Lmod1_peaks <- mincFindPeaks(Model1, column = 'tvalue-poly(Age, 2)1:GroupAlc', direction = "both", threshold = 3.8, minDistance = 4.5) 
#Pmod1_peaks <- mincFindPeaks(Model1, column = 'tvalue-poly(Age, 2)2:GroupAlc', direction = "both", threshold = 3, minDistance = 6) # non significant after correction

Lmod2H_peaks <- mincFindPeaks(Model2, column = 'tvalue-poly(Age, 2)1:IntakeHigh', direction = "both", threshold = 3.1, minDistance = 8) 
Lmod2L_peaks <- mincFindPeaks(Model2, column = 'tvalue-poly(Age, 2)1:IntakeLow', direction = "both", threshold = 4.48, minDistance = 4) 
Pmod2H_peaks <- mincFindPeaks(Model2, column = 'tvalue-poly(Age, 2)2:IntakeHigh', direction = "both", threshold = 4.39, minDistance = 8) 
#Pmod2L_peaks <- mincFindPeaks(Model2, column = 'tvalue-poly(Age, 2)2:IntakeLow', direction = "both", threshold = 3, minDistance = 8) # non significant

Lmod3_peaks <- mincFindPeaks(Model3, column = 'tvalue-poly(Age, 2)1:GroupAlc:Sexmale', direction = "both", threshold = 3.71, minDistance = 4) 
Pmod3_peaks <- mincFindPeaks(Model3, column = 'tvalue-poly(Age, 2)2:GroupAlc:Sexmale', direction = "both", threshold = 3.74, minDistance = 4) 

Lmod4H_peaks <- mincFindPeaks(Model4, column = 'tvalue-poly(Age, 2)1:IntakeHigh:Sexmale', direction = "both", threshold = 4.06, minDistance = 8) 
Lmod4L_peaks <- mincFindPeaks(Model4, column = 'tvalue-poly(Age, 2)1:IntakeLow:Sexmale', direction = "both", threshold = 3.51, minDistance = 8) 
Pmod4H_peaks <- mincFindPeaks(Model4, column = 'tvalue-poly(Age, 2)2:IntakeHigh:Sexmale', direction = "both", threshold = 3.98, minDistance = 8)
Pmod4L_peaks <- mincFindPeaks(Model4, column = 'tvalue-poly(Age, 2)2:IntakeLow:Sexmale', direction = "both", threshold = 3.73, minDistance = 4) 

Lmod5_peaks <- mincFindPeaks(Model5, column = 'tvalue-poly(Age, 2)1:IntakeHigh', direction = "both", threshold = 3.34, minDistance = 10) 
Pmod5_peaks <- mincFindPeaks(Model5, column = 'tvalue-poly(Age, 2)2:IntakeHigh', direction = "both", threshold = 3.61, minDistance = 8)

Lmod6_peaks <- mincFindPeaks(Model6, column = 'tvalue-poly(Age, 2)1:IntakeHigh:Sexmale', direction = "both", threshold = 3.83, minDistance = 7) 
Pmod6_peaks <- mincFindPeaks(Model6, column = 'tvalue-poly(Age, 2)2:IntakeHigh:Sexmale', direction = "both", threshold = 3.77, minDistance = 7) 

mod_peaks <- list(Lmod1_peaks,Lmod2H_peaks,Lmod2L_peaks,Pmod2H_peaks,Lmod3_peaks,Pmod3_peaks,Lmod4H_peaks,Lmod4L_peaks,Pmod4H_peaks,Pmod4L_peaks,Lmod5_peaks,Pmod5_peaks,Lmod6_peaks,Pmod6_peaks) %>% 
  map(~ mincLabelPeaks(.x,  atlas, defs="data/SIGMA_InVivo_Anatomical_Brain_Atlas_ListOfStructures_dav.csv") ) %>% 
  set_names(c("Lmod1_peaks","Lmod2H_peaks","Lmod2L_peaks","Pmod2H_peaks","Lmod3_peaks","Pmod3_peaks","Lmod4H_peaks","Lmod4L_peaks","Pmod4H_peaks","Pmod4L_peaks",
	      "Lmod5_peaks","Pmod5_peaks","Lmod6_peaks","Pmod6_peaks"))

dir.create("Peaks")

mod_peaks %>% iwalk(~write_csv(.x, paste0(getwd(),"/Peaks/",.y, ".csv")))

# Trayectories -------------------------------------------------------------

mod_peaks <- mod_peaks %>% map(~ .x %>% mutate(label_clean = label %>% make_clean_names()) )

ROIs <- names(mod_peaks) %>% map(~ mod_peaks[[.x]] %>% add_column(Model = rep(.x, nrow(mod_peaks[[.x]]) )) ) %>% reduce(rbind) %>% mutate(ROI_model = str_c(label_clean,'_',Model) )

Jdata_jacobians <- Jdata %>% select(RID,Session,Subject,Group,Intake,Age,Sex,Batch)

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
                                                Group == "Ctrl" & Sex == "male" ~ "Ctrl M"),
			  Intake_sex = case_when(Intake == "High" & Sex == "female" ~ "High F",
                                                 Intake == "High" & Sex == "male" ~ "High Me",
						 Intake == "Low" & Sex == "female" ~ "Low F",
                                                 Intake == "Low" & Sex == "male" ~ "Low M",
                                                 Intake == "Ctrl" & Sex == "female" ~ "Ctrl F",
                                                 Intake == "Ctrl" & Sex == "male" ~ "Ctrl M"), .before = 8) %>%
			  mutate(Group_sex = factor(Group_sex, levels = c("Ctrl F","Ctrl M","Alc F","Alc M") ),
			  	 Intake_sex = factor(Intake_sex, levels = c("Ctrl F","Ctrl M","Low F","Low M","High F","High M")) )

# Plotting ----------------------------------------------------------------

theme_settings <- theme(text = element_text(size=20),
                        axis.text.x = element_text(size=15))

# colors
pal_group <- c("#737373","#B22222","#217175")
pal_sex <- c(alpha("#8d289f",1),alpha("#28799f",1))
pal_groupSex <- c("#737373","#424242","#008B8B","#8B475D")
pal_IntakeSex <- c("#737373","#424242","#7DA2E8","#00688B","#9780F4","#6C3F81")
                        
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
		  aes(y=fit, ymin=lower, ymax=upper, fill = Group), alpha=0.3) +
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

# Intake

dir.create("Trayectories/Linear/Intake", recursive = TRUE)

plots_ROIs_Intake_L <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_Intake_L[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Intake",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Intake", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Intake", "Age"),
                                          lmer(get(ROI) ~ Age*Intake + 
                                                 Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Intake), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Intake", "Age"),
                                            lmer(get(ROI) ~ Age*Intake + 
                                                   Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                            xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                                 max(Jdata_jacobians$Age),1)))), aes(y=fit, ymin=lower,
                                                                                                     ymax=upper, fill = Intake), alpha=0.3) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_Intake_L) %>% map(~ ggsave(filename = paste0("Trayectories/Linear/Intake/",plots_ROIs_Intake_L[.x] %>% names(),".png"), plot = plots_ROIs_Intake_L[[.x]],dpi = 300, width = 5.5, height = 4.5) )


dir.create("Trayectories/Linear/Intake_sex", recursive = TRUE)

plots_ROIs_Intake_sex_L <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_Intake_sex_L[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Intake_sex",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Intake_sex", linewidth =2,
                  palette = pal_IntakeSex,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Intake_sex", "Age"),
                                          lmer(get(ROI) ~ Age*Intake_sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Intake_sex), size=2) +
	theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_Intake_sex_L) %>% map(~ ggsave(filename = paste0("Trayectories/Linear/Intake_sex/",plots_ROIs_Intake_sex_L[.x] %>% names(),".png"), plot = plots_ROIs_Intake_sex_L[[.x]],dpi = 300, width = 5.5, height = 4.5) )


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

# Intake

dir.create("Trayectories/Polynomial/Intake", recursive = TRUE)

plots_ROIs_Intake <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_Intake[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Intake",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Intake", linewidth =2,
                  palette = pal_group,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Intake", "Age"),
                                          lmer(get(ROI) ~ poly(Age,2)*Intake + 
                                                 Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Intake), size=2) +
        geom_ribbon(data = as_tibble(Effect(c("Intake", "Age"),
                                            lmer(get(ROI) ~ poly(Age,2)*Intake + 
                                                   Sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                            xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                                 max(Jdata_jacobians$Age),1)))), aes(y=fit, ymin=lower,
                                                                                                     ymax=upper, fill = Intake), alpha=0.3) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_Intake) %>% map(~ ggsave(filename = paste0("Trayectories/Polynomial/Intake/",plots_ROIs_Intake[.x] %>% names(),".png"), plot = plots_ROIs_Intake[[.x]],dpi = 300, width = 5.5, height = 4.5) )


dir.create("Trayectories/Polynomial/Intake_sex", recursive = TRUE)

plots_ROIs_Intake_sex <- NULL
for (i in 1:length(nROI)) {
  ROI <- nROI[i]
  plots_ROIs_Intake_sex[[ROI]] <- ggscatter(Jdata_jacobians,
                  x = "Age", y = ROI, group = "Intake_sex",
                  #add = c("loess"), conf.int = TRUE,
                  color = "Intake_sex", linewidth =2,
                  palette = pal_IntakeSex,
                  add.params = list(size = 2, alpha = 0.5),
                  plot_type ="b",
                  title = gsub('_left','',ROIs$label[i]),
                  xlab = "Age (PND)",
                  ylab = "Local volume",
                  font.x = c(16,"bold"),
                  font.y = c(16,"bold"),
                  font.tickslab = c(14,"bold")) +
        geom_line(data = as_tibble(Effect(c("Intake_sex", "Age"),
                                          lmer(get(ROI) ~ poly(Age,2)*Intake_sex + Batch + (1 |RID), data = Jdata_jacobians), 
                                          xlevels=list(Age=seq(min(Jdata_jacobians$Age),
                                                               max(Jdata_jacobians$Age),1)))), 
                  aes(y=fit, color = Intake_sex), size=2) +
	theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,size = 16, face = "bold")) + 
        theme_settings 
}

1:length(plots_ROIs_Intake_sex) %>% map(~ ggsave(filename = paste0("Trayectories/Polynomial/Intake_sex/",plots_ROIs_Intake_sex[.x] %>% names(),".png"), plot = plots_ROIs_Intake_sex[[.x]],dpi = 300, width = 5.5, height = 4.5) )


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
                      Model2 = case_when(Model == "Lmod1_peaks" ~ "Model 1",
                                        Model == "Lmod2H_peaks" ~ "Model 2",
                                        Model == "Lmod2L_peaks" ~ "Model 2",
                                        Model == "Lmod3_peaks" ~ "Model 3",
                                        Model == "Lmod4H_peaks" ~ "Model 4",
                                        Model == "Lmod4L_peaks" ~ "Model 4",
                                        Model == "Lmod5_peaks" ~ "Model 5",
                                        Model == "Lmod6_peaks" ~ "Model 6",
					Model == "Pmod1_peaks" ~ "Model 1",
                                        Model == "Pmod2H_peaks" ~ "Model 2",
                                        Model == "Pmod2L_peaks" ~ "Model 2",
                                        Model == "Pmod3_peaks" ~ "Model 3",
                                        Model == "Pmod4H_peaks" ~ "Model 4",
                                        Model == "Pmod4L_peaks" ~ "Model 4",
                                        Model == "Pmod5_peaks" ~ "Model 5",
                                        Model == "Pmod6_peaks" ~ "Model 6"),
                      Contrast = case_when(Model == "Lmod1_peaks" ~ "EtOH > Ctrl",
                                        Model == "Lmod2H_peaks" ~ "High EtOH > Ctrl",
                                        Model == "Lmod2L_peaks" ~ "Low EtOH > Ctrl",
                                        Model == "Lmod3_peaks" ~ "EtOH > Ctrl",
                                        Model == "Lmod4H_peaks" ~ "High EtOH > Ctrl",
                                        Model == "Lmod4L_peaks" ~ "Low EtOH > Ctrl",
                                        Model == "Lmod5_peaks" ~ "High EtOH > Low",
                                        Model == "Lmod6_peaks" ~ "High EtOH > Low",
					Model == "Pmod1_peaks" ~ "EtOH > Ctrl",
                                        Model == "Pmod2H_peaks" ~ "High EtOH > Ctrl",
                                        Model == "Pmod2L_peaks" ~ "Low EtOH > Ctrl",
                                        Model == "Pmod3_peaks" ~ "EtOH > Ctrl",
                                        Model == "Pmod4H_peaks" ~ "High EtOH > Ctrl",
                                        Model == "Pmod4L_peaks" ~ "Low EtOH > Ctrl",
                                        Model == "Pmod5_peaks" ~ "High EtOH > Low",
                                        Model == "Pmod6_peaks" ~ "High EtOH > Low"),
 		      Type = case_when(Model == "Lmod1_peaks" ~ "Linear",
                                        Model == "Lmod2H_peaks" ~ "Linear",
                                        Model == "Lmod2L_peaks" ~ "Linear",
                                        Model == "Lmod3_peaks" ~ "Linear",
                                        Model == "Lmod4H_peaks" ~ "Linear",
                                        Model == "Lmod4L_peaks" ~ "Linear",
                                        Model == "Lmod5_peaks" ~ "Linear",
                                        Model == "Lmod6_peaks" ~ "Linear",
					Model == "Pmod1_peaks" ~ "Poly",
                                        Model == "Pmod2H_peaks" ~ "Poly",
                                        Model == "Pmod2L_peaks" ~ "Poly",
                                        Model == "Pmod3_peaks" ~ "Poly",
                                        Model == "Pmod4H_peaks" ~ "Poly",
                                        Model == "Pmod4L_peaks" ~ "Poly",
                                        Model == "Pmod5_peaks" ~ "Poly",
                                        Model == "Pmod6_peaks" ~ "Poly") ) %>% 
  select(Model2,Type,Contrast,ROIs, Hemisphere, Coordinates, Effect, "t-value")

write_csv(ROIs,"Trayectories/ROIs.csv")
write_csv(ROI_table,"Trayectories/ROI_table.csv")

save(atlas,mod_peaks,theme_settings,plots_ROIs_L,plots_ROIs_sex_L,Jdata_jacobians,ROIs,ROI_table,nROI, file = "Trayectories_data1.RData")
save(plots_ROIs_Intake_L,plots_ROIs_Intake_sex_L, file = "Trayectories_data2.RData")
save(plots_ROIs,plots_ROIs_sex, file = "Trayectories_data3.RData")
save(plots_ROIs_Intake,plots_ROIs_Intake_sex, file = "Trayectories_data4.RData")
