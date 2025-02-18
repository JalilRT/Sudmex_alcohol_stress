library(tidyverse)

path="/scratch/m/mchakrav/aguilars/MethAgg_project/batch1/MRI"
setwd(path)

# Create directories

full_path=list.files(path = paste0(path,"/jacobians/output/secondlevel/resampled-dbm/jacobian/relative/smooth/"), 
    pattern="_1mm.mnc", recursive=T, full.names=T) %>% as_tibble()

detect_files <- full_path %>% mutate(file = value %>% basename(),
  RID = file %>% str_split(pattern="_") %>% map(~ .x %>% .[1]) %>% reduce(rbind),
  Session = file %>% str_split(pattern="_") %>% map(~ .x %>% .[2]) %>% reduce(rbind), .before=1 ) %>%
  rename("Subject" =	"value") %>% select(-(file)) %>% add_column(Sex = "female", .after = 2)

RID <- read_table(paste0(path, "/code/RID.csv")) %>% 
    mutate(RID = case_when(
            as.numeric(RID) <= 10 ~ paste0("sub-00", RID),
            as.numeric(RID) <= 99 ~ paste0("sub-0", RID),
            TRUE ~ paste0("sub-0", RID)),
        Session = case_when(
            Sesion == "T1" ~ paste0("ses-", Sesion),
            Sesion == "T2" ~ paste0("ses-", Sesion) )) %>% select(- c(Sesion))

fc_input <- left_join(detect_files,RID,by=c("RID","Session")) %>% 
    mutate_all(as.character) %>% 
    mutate_all(~ ifelse(is.na(.), "", .)) %>% 
    mutate_all(~ gsub(",", "", .))

write_csv(fc_input, paste0(path, "/code/DBM_dataset.csv"), col_names = TRUE)

