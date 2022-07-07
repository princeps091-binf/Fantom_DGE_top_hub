library(edgeR)
library(tidyverse)
library(viridis)
library(fgsea)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(out_tbl)
}


#-----------------------------------------

count_data_tbl<-tbl_in_fn("./data/enh_CAGE_count_tbl.Rda")
count_data_tbl<-count_data_tbl %>% 
  filter(if_any(where(is.numeric), ~ .x > 1))


group<-as.factor(rep(c("HMEC","MCF7","MDA"),c(3,4,1)))
design <- model.matrix(~group)

y <- DGEList(counts=count_data_tbl[,-1],genes=count_data_tbl[,1])
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

qlf.MCF7vsHMEC <- glmQLFTest(fit, coef=2)
qlf.MDAvsHMEC <- glmQLFTest(fit, coef=3)


topTags(qlf.MCF7vsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(logFC))+
  geom_density()

topTags(qlf.MCF7vsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(PValue))+
  geom_histogram(bins=100)

topTags(qlf.MCF7vsHMEC,n = Inf)[[1]] %>% 
  as_tibble %>% 
  ggplot(.,aes(logFC,-log10(FDR)))+
  geom_point(size=0.1,alpha=0.2)
