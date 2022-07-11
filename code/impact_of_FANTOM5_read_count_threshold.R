library(DESeq2)
library(tidyverse)
library(viridis)
library(kableExtra)
library(fgsea)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(out_tbl)
}


#-----------------------------------------

#load the count datasets
##HMEC
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_hmec_count_tbl.Rda")
cage_count_hmec<-cage_count_hmec%>%dplyr::slice(-1)
##MDA
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_mda_count_tbl.Rda")
cage_count_MDA<-cage_count_MDA%>%dplyr::slice(-1)
##MCF7
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_mcf7_count_tbl.Rda")
cage_count_mcf7<-cage_count_mcf7%>%dplyr::slice(-1)

colnames(cage_count_hmec)<-c("ID","HMEC.1","HMEC.2","HMEC.3","HMEC.m")
colnames(cage_count_MDA)<-c("ID","MDA.1","MDA.m")
colnames(cage_count_mcf7)<-c("ID","MCF7.1","MCF7.2","MCF7.3","MCF7.4","MCF7.m")


count_data_tbl<-cage_count_hmec%>%dplyr::select(-contains("m",ignore.case = F))%>%
  full_join(.,cage_count_MDA%>%dplyr::select(-contains("m",ignore.case = F)))%>%
  full_join(.,cage_count_mcf7%>%dplyr::select(-contains("m",ignore.case = F)))
count_data_tbl<-count_data_tbl %>% replace(is.na(.), 0)

fantom_peak_count<-apply(count_data_tbl[,-1],1,function(x){
  any(x>1)})

current_peak_count<-nrow(count_data_tbl %>% 
  mutate(min.count=fantom_peak_count))
fantom_peak_count<-nrow(count_data_tbl %>% 
       mutate(min.count=fantom_peak_count) %>% 
       filter(fantom_peak_count))

enh_count_data_tbl<-tbl_in_fn("./data/enh_CAGE_count_tbl.Rda")
fantom_enh_peak_count<-apply(enh_count_data_tbl[,-1],1,function(x){
  any(x>1)})

current_enh_count<-nrow(enh_count_data_tbl %>% 
                           mutate(min.count=fantom_enh_peak_count))
fantom_enh_count<-nrow(enh_count_data_tbl %>% 
                          mutate(min.count=fantom_enh_peak_count) %>% 
                          filter(fantom_enh_peak_count))


tibble(peak.count=c(current_peak_count,fantom_peak_count,current_enh_count,fantom_enh_count),
       set=rep(c('Current (>0 read)','FANTOM5 (>1 read)'),2),
       enh.tss=rep(c("TSS","Enhancer"),each=2)) %>% 
  ggplot(.,aes(set,peak.count))+
  geom_bar(stat='identity')+
  theme_minimal()+
  facet_grid(enh.tss~.,scales="free")
ggsave("~/Documents/multires_bhicect/weeklies/weekly61/img/FANTOM5_criteria_effect.svg")
#-------------------
tss_dge_tbl<-tbl_in_fn("./data/tss_2tag_HMEC_MCF7_MDA_edgeR.Rda")
enh_dge_tbl<-tbl_in_fn("./data/enh_2tag_HMEC_MCF7_MDA_edgeR.Rda")

tss_dge_tbl %>% 
  ggplot(.,aes(mcf7.pval))+
  geom_histogram(bins=100)
ggsave("~/Documents/multires_bhicect/weeklies/weekly61/img/FANTOM5_criteria_tss_pval_effect_edgeR.svg")

enh_dge_tbl %>% 
  ggplot(.,aes(mcf7.pval))+
  geom_histogram(bins=100)
ggsave("~/Documents/multires_bhicect/weeklies/weekly61/img/FANTOM5_criteria_enh_pval_effect_edgeR.svg")
