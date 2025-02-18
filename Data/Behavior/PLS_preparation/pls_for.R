setwd("/Users/jalil/phd/PhD/Psilantro/Sudmex_alcohol_stress/Data/Behavior/PLS_preparation/")

library(tidyverse)
library(magrittr)

data <- read_csv("Behavior_metrics4pls.csv") %>% 
  select(-Batch) %>% 
  #replace -0 with - in RID column
  mutate(RID = str_replace_all(RID,"-0","-"))

files_names <- list.files("/Users/jalil/phd/PhD/Psilantro/Sudmex_alcohol_stress/Data/MRI/DBM/jacobians/") %>% 
  #keep only those who has "T3" in the name
  keep(~ grepl("T3", .)) %>%
  keep(~ grepl("nii.gz", .)) 

files_names %>%
  str_split(pattern="_") %>% map(~ .x[[1]]) %>% reduce(rbind) %>% as_tibble() %>% 
  add_column(relative_jacobian = files_names) %>% 
  set_colnames(c("RID","relative_jacobian")) %>% right_join(data,by = "RID") %>% 
  mutate(relative_jacobian = paste0("/data/chamal/projects/jalilr/pls/data/relative_jacobians/nii/t3_pls/",relative_jacobian)) %>% 
  write_csv("Behavior_4pls_path.csv")

paste(colnames(data) %>% c() %>% noquote() %>% .[-c(1,2,4)],collapse = "','")

      