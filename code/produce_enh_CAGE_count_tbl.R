library(tidyverse)
library(vroom)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------

enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt")

MDA<-"CNhs10736"
MCF7<-c("CNhs11943","CNhs12564","CNhs12475","CNhs12703")
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')

enh_count_tbl<-enh_tbl %>% 
  dplyr::select(Id,c(HMEC,MCF7,MDA)) %>% 
  dplyr::rename("HMEC.1"="CNhs11077","HMEC.2"="CNhs11382","HMEC.3"="CNhs12032") %>% 
  dplyr::rename("MCF7.1"="CNhs11943","MCF7.2"="CNhs12564","MCF7.3"="CNhs12475","MCF7.4"="CNhs12703") %>% 
  dplyr::rename("MDA.1"="CNhs10736") %>% 
  filter(if_any(where(is.numeric), ~ .x > 0))
  
save(enh_count_tbl,file = "./data/enh_CAGE_count_tbl.Rda")
