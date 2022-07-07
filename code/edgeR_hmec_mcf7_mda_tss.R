library(edgeR)
library(tidyverse)
library(viridis)
library(fgsea)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
#load the count datasets
##HMEC
base::load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_hmec_count_tbl.Rda")
cage_count_hmec<-cage_count_hmec%>%dplyr::slice(-1)
##MDA
base::load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_mda_count_tbl.Rda")
cage_count_MDA<-cage_count_MDA%>%dplyr::slice(-1)
##MCF7
base::load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_mcf7_count_tbl.Rda")
cage_count_mcf7<-cage_count_mcf7%>%dplyr::slice(-1)

colnames(cage_count_hmec)<-c("ID","HMEC.1","HMEC.2","HMEC.3","HMEC.m")
colnames(cage_count_MDA)<-c("ID","MDA.1","MDA.m")
colnames(cage_count_mcf7)<-c("ID","MCF7.1","MCF7.2","MCF7.3","MCF7.4","MCF7.m")


count_data_tbl<-cage_count_hmec%>%dplyr::select(-contains("m",ignore.case = F))%>%
  full_join(.,cage_count_MDA%>%dplyr::select(-contains("m",ignore.case = F)))%>%
  full_join(.,cage_count_mcf7%>%dplyr::select(-contains("m",ignore.case = F)))
count_data_tbl<-count_data_tbl %>% replace(is.na(.), 0)
count_data_tbl<-count_data_tbl %>% 
  filter(if_any(where(is.numeric), ~ .x > 1))

group<-as.factor(rep(c("HMEC","MDA","MCF7"),c(3,1,4)))
design <- model.matrix(~group)

y <- DGEList(counts=count_data_tbl[,-1],genes=count_data_tbl[,1])
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

qlf.MCF7vsHMEC <- glmQLFTest(fit, coef=2)
qlf.MDAvsHMEC <- glmQLFTest(fit, coef=3)


topTags(qlf.MDAvsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(logFC))+
  geom_density()

topTags(qlf.MDAvsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(PValue))+
  geom_histogram(bins=100)

topTags(qlf.MDAvsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(logFC,-log10(FDR)))+
  geom_point(size=0.1)
