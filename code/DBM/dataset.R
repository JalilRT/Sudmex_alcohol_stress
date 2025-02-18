setwd("/home/jalil/Documents/borrar/rabies/stress/")

pacman::p_load(tidyverse,readxl,magrittr)

rq <- read_excel("tividavo.xlsx")
qr <- read_csv("RID.csv")

rsq <- rq %>% mutate(RID=subjects %>% basename() %>% 
                str_split(pattern = "_") %>% map(~ .x %>% .[[2]]) %>% reduce(rbind),
              Ses=subjects %>% basename() %>% 
                str_split(pattern = "_") %>% map(~ .x %>% .[[3]]) %>% reduce(rbind),.before = 1) %>% 
  set_colnames(c("RID","Ses","subjects"))

rsq1 <- rsq %>% filter(Ses=="ses-T1") %>% left_join(qr,by="RID")
rsq2 <- rsq %>% filter(Ses=="ses-T2") %>% left_join(qr,by="RID")
rsq3 <- rsq %>% filter(Ses=="ses-T3") %>% left_join(qr,by="RID")
rsq5 <- rsq %>% filter(Ses=="ses-T5") %>% left_join(qr,by="RID")

rsq <- rbind(rsq1,rsq2,rsq3,rsq5)
write.csv(rsq,"cuca.csv")

