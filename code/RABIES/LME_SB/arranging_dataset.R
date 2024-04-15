library(tidyverse)

path="/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri"
setwd(path)

fc_csv <- read_csv("NBR/fc_dataset.csv")

# Create directories
c(1,2,3,5) %>% map(~ dir.create(paste0("Seed_based/ses-T",.x,"/Alc"), recursive=T) )
c(1,2,3,5) %>% map(~ dir.create(paste0("Seed_based/ses-T",.x,"/Ctrl"), recursive=T) )

# list rois
rois=c("left_granule_cell_level_of_the_cerebellum_left_2","left_granule_cell_level_of_the_cerebellum_left_Pmod4L","left_molecular_layer_of_the_cerebellum",
	"left_pre_limbic_system_left","left_presubiculum_left_Pmod6","right_cornu_ammonis_2","right_granule_cell_level_of_the_cerebellum",
	"right_hypothalamic_region","right_olfactory_bulb")

for(roi in rois){
 # Arranging
 full_path=list.files(path = paste0(path,"/Seed_based"), pattern=roi, recursive=T, full.names=T) %>% as_tibble()

 detect_files <- full_path %>% mutate(file = value %>% basename(),
   RID = file %>% str_split(pattern="_") %>% map(~ .x %>% .[1]) %>% reduce(rbind),
   Ses = file %>% str_split(pattern="_") %>% map(~ .x %>% .[2]) %>% reduce(rbind), .before=1 ) %>%
   rename("InputFile" =	"value") %>% select(-(file))

 fc_input <- left_join(fc_csv,detect_files,by=c("RID","Ses")) %>% rename("Subj" = "RID")

 write.table(x = fc_input, file = paste0("code/LME_SB/fc_dataset_",roi,".txt"), sep = " ", row.names=FALSE, quote=FALSE)
}
