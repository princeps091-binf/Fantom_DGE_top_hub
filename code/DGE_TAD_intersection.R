library(GenomicRanges)
library(tidyverse)
library(furrr)
library(DESeq2)
#---------------------------------------------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#-----------------------------------------
TAD_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_H1_TAD_pval_tbl.Rda"
dge_file<-"./data/HMEC_MCF7_MDA_DSeq2.Rda"

#-----------------------------------------

TAD_tbl<-get_obj_in_fn(TAD_file)%>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01)
TAD_GRange<-IRanges::reduce(do.call("c",unlist(TAD_tbl$GRange)))


dge_tbl<-get_obj_in_fn(dge_file)

res_mcf7 <- results( dge_tbl ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( dge_tbl ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(res_mcf7_tbl,res_mda_tbl,res_mda,res_mcf7)
rm(dge_tbl)
dge_coord_tbl<-res_dge_tbl %>%
  dplyr::select(ID) %>% 
  mutate(chr=str_split_fixed(ID,":|\\.\\.|,",4)[,1],
         start=as.numeric(str_split_fixed(ID,":|\\.\\.|,",4)[,2]),
         end=as.numeric(str_split_fixed(ID,":|\\.\\.|,",4)[,3]))
dge_Grange<-   GRanges(seqnames=dge_coord_tbl$chr,
                       ranges = IRanges(start=dge_coord_tbl$start,
                                        end=dge_coord_tbl$end
                       ))
mcols(dge_Grange)<-res_dge_tbl

ok_dge_GRange<-dge_Grange[unique(c(which(!(is.na(mcols(dge_Grange)$mcf7.padj))),which(!(is.na(mcols(dge_Grange)$mda.padj)))))]

in_peak<-mcols(ok_dge_GRange)$ID[unique(queryHits(findOverlaps(ok_dge_GRange,TAD_GRange)))]

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_peak,"in","out")) %>% 
  ggplot(.,aes(mda.lfc,-log10(mda.padj)))+
  geom_point(alpha=0.1)+
  facet_grid(hub.io~.,scales="free")

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_peak,"in","out")) %>% 
  ggplot(.,aes(mcf7.lfc,-log10(mcf7.padj)))+
  geom_point(alpha=0.1)+
  facet_grid(hub.io~.,scales="free")

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_peak,"in","out")) %>% 
  ggplot(.,aes(mda.lfc,color=hub.io))+
  geom_density()


