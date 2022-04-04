#Compute DGE 
library(DESeq2)
library(tidyverse)
library(viridis)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

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

col_data_tbl<-tibble(cell.line=as.factor(rep(c("HMEC","MDA","MCF7"),c(3,1,4))))

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count_data_tbl[,-1],
  colData = col_data_tbl,
  design = ~ cell.line)
rownames(ddsFullCountTable)<-count_data_tbl$ID

dds <- DESeq(ddsFullCountTable)
save(dds,file = "./data/HMEC_MCF7_MDA_DSeq2.Rda")
resultsNames(dds)

res_mcf7 <- results( dds ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( dds ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(col_data_tbl,count_data_tbl,ddsFullCountTable)
rm(cage_count_hmec,cage_count_mcf7,cage_count_MDA)
rm(res_mcf7_tbl,res_mda_tbl,res_mda,res_mcf7)
rm(dds)

res_dge_tbl %>% arrange(desc(mcf7.lfc)) %>% 
  mutate(idx=1:n(),signif=ifelse(is.na(mcf7.padj),"out",ifelse(mcf7.padj<0.01,"FDR < 0.01","FDR > 0.01"))) %>% 
#  filter(sign(mcf7.lfc+mcf7.lfc.se)*sign(mcf7.lfc-mcf7.lfc.se)==1) %>% 
  ggplot(aes(x=mcf7.lfc-mcf7.lfc.se,xend=mcf7.lfc+mcf7.lfc.se,y=idx,yend=idx,color=signif))+
  geom_segment(alpha=1,size=0.05)+
  geom_vline(xintercept = 0,linetype=2)+facet_wrap(signif~.,scales="free")+
  theme(legend.position = "none")

res_dge_tbl %>% arrange(desc(mda.lfc)) %>% 
  mutate(idx=1:n(),signif=ifelse(is.na(mda.padj),"out",ifelse(mda.padj<0.01,"FDR < 0.01","FDR > 0.01"))) %>% 
#  filter(sign(mda.lfc+mda.lfc.se)*sign(mda.lfc-mda.lfc.se)==1) %>% 
  ggplot(aes(x=mda.lfc-mda.lfc.se,xend=mda.lfc+mda.lfc.se,y=idx,yend=idx,color=signif))+
  geom_segment(alpha=0.5,size=0.05)+
  geom_vline(xintercept = 0,linetype=2)+facet_wrap(signif~.)+
  theme(legend.position = "none")

res_dge_tbl %>% 
  mutate(signif=ifelse(is.na(mda.padj) | is.na(mcf7.padj),"out","in")) %>% 
  ggplot(.,aes(mcf7.lfc,mda.lfc))+
  geom_point(alpha=0.01)+
  facet_grid(signif~.,scales='free')

do.call(bind_rows,map(seq(0,1,length.out=501),function(x){
  return(res_dge_tbl %>% 
    filter(mcf7.padj<=x | mda.padj<=x) %>% 
    mutate(dge.ok=as.character(ifelse(sign(mcf7.lfc)*sign(mda.lfc)>0,"pro","contra"))) %>% 
    group_by(dge.ok) %>% 
    summarise(n=n()) %>% 
    mutate(FDR=x))
  
})) %>% 
  ggplot(.,aes(FDR,n,fill=dge.ok))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")
