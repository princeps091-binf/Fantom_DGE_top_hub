library(DESeq2)
library(tidyverse)
library(vroom)
library(furrr)

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

count_data_tbl<-tbl_in_fn("./data/enh_CAGE_count_tbl.Rda")
count_data_tbl<-count_data_tbl %>% 
  filter(if_any(where(is.numeric), ~ .x > 1))

col_data_tbl<-tibble(cell.line=as.factor(rep(c("HMEC","MCF7","MDA"),c(3,4,1))))

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count_data_tbl[,-1],
  colData = col_data_tbl,
  design = ~ cell.line)
rownames(ddsFullCountTable)<-count_data_tbl$Id

dds <- DESeq(ddsFullCountTable)
resultsNames(dds)
res_mcf7 <- results( dds ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( dds ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(col_data_tbl,count_data_tbl,ddsFullCountTable)

save(res_dge_tbl,file="./data/enh_HMEC_MCF7_MDA_DSeq2.Rda")
