library(tidyverse)
library(GenomicRanges)
library(DESeq2)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(out_tbl)
}

#-------------------------------------------------------------------------------------------------------
compound_hub_5kb_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
TAD_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_TAD_pval_tbl.Rda"

dge_file<-"./data/HMEC_MCF7_MDA_DSeq2.Rda"

compound_hub_5kb_tbl<-tbl_in_fn(compound_hub_5kb_file)


top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))

top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  tmp_tbl<-top_compound_hub_5kb_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
    mutate(bins=map(bins,as.numeric)) 
  
})) %>% 
  dplyr::select(chr,parent.hub,res,bins)

top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(GRange=pmap(list(chr,bins,res),function(chr,bins,res){
    inter_cl_Grange<-   GRanges(seqnames=chr,
                                ranges = IRanges(start=as.numeric(bins),
                                                 end=as.numeric(bins) + res_num[res]-1
                                ))
    inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
    return(inter_cl_Grange)
    
  }))
cl_GRange<-IRanges::reduce(do.call("c",top_compound_hub_5kb_tbl$GRange))
#-------------------------------------
TAD_tbl<-tbl_in_fn(TAD_file)%>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01)
TAD_GRange<-IRanges::reduce(do.call("c",unlist(TAD_tbl$GRange)))
#-------------------------------------

dge_tbl<-tbl_in_fn(dge_file)

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

in_cl_peak<-mcols(ok_dge_GRange)$ID[unique(queryHits(findOverlaps(ok_dge_GRange,cl_GRange)))]

in_tad_peak<-mcols(ok_dge_GRange)$ID[unique(queryHits(findOverlaps(ok_dge_GRange,TAD_GRange)))]

library(UpSetR)
up_l<-list(TAD=in_tad_peak,hub=in_cl_peak)
UpSetR::upset(fromList(up_l))

tad_ex_peak<-fromList(up_l) %>%
  mutate(peak=unique(unlist(up_l))) %>% 
  filter(TAD==1 & hub ==0) %>% 
  dplyr::select(peak) %>% unlist

peak_ex_peak<-fromList(up_l) %>%
  mutate(peak=unique(unlist(up_l))) %>% 
  filter(TAD==0 & hub ==1) %>% 
  dplyr::select(peak) %>% unlist

inter_peak<-fromList(up_l) %>%
  mutate(peak=unique(unlist(up_l))) %>% 
  filter(TAD==1 & hub ==1) %>% 
  dplyr::select(peak) %>% unlist

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_tad_peak,"in","out")) %>% 
  ggplot(.,aes(mda.lfc,-log10(mda.padj)))+
  geom_point(alpha=0.1)+
  facet_grid(hub.io~.,scales="free")

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% inter_peak,"inter",ifelse(ID %in% tad_ex_peak,"tad",ifelse(ID %in% peak_ex_peak,"hub","out")))) %>% 
  ggplot(.,aes(mda.lfc,-log10(mda.padj)))+
  geom_point(alpha=0.01)+
  facet_grid(hub.io~.,scales="free")

res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% inter_peak,"inter",ifelse(ID %in% tad_ex_peak,"tad",ifelse(ID %in% peak_ex_peak,"hub","out")))) %>% 
  filter(hub.io!="out") %>% 
  ggplot(.,aes(mda.lfc,color=hub.io))+
  geom_density()
