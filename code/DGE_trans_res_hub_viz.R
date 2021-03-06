library(tidyverse)
library(GenomicRanges)
library(DESeq2)
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

#-------------------------------------------------------------------------------------------------------
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
dge_file<-"./data/tss_2tag_HMEC_MCF7_MDA_edgeR.Rda"
#-------------------------------------------------------------------------------------------------------
hub_tbl<-tbl_in_fn(hub_file) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

hub_tbl<-do.call(bind_rows,map(unique(hub_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  tmp_tbl<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,as.numeric)) 
  
})) %>% 
  dplyr::select(chr,node,res,bins)

hub_tbl<-hub_tbl %>% 
  mutate(GRange=pmap(list(chr,bins,res),function(chr,bins,res){
    inter_cl_Grange<-   GRanges(seqnames=chr,
                                ranges = IRanges(start=as.numeric(bins),
                                                 end=as.numeric(bins) + res_num[res]-1
                                ))
    inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
    return(inter_cl_Grange)
    
  }))
cl_GRange<-IRanges::reduce(do.call("c",hub_tbl$GRange))
tmp_l<-hub_tbl %>% 
  #  filter(res%in% c("1Mb","500kb","100kb")) %>% 
  #  filter(res%in% c("10kb","50kb","5kb")) %>% 
  dplyr::select(GRange) %>% as.list
cl_GRange<-IRanges::reduce(do.call("c",tmp_l$GRange))
#-------------------------------------------------------------------------------------------------------
dge_tbl<-tbl_in_fn(dge_file)

res_mcf7 <- results( dge_tbl ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( dge_tbl ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(res_mcf7_tbl,res_mda_tbl,res_mda,res_mcf7)
rm(dge_tbl)

res_dge_tbl<-tbl_in_fn(dge_file)

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


gg_tmp<-res_dge_tbl %>% 
#  filter(ID %in% mcols(ok_dge_GRange)$ID) %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  ggplot(.,aes(mcf7.lfc,-log10(mcf7.padj),color=hub.io))+
  geom_point(size=0.01)+
  facet_grid(hub.io~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly59/img/volcano_MCF7.png",gg_tmp)
gg_tmp<-res_dge_tbl %>% 
#  filter(ID %in% mcols(ok_dge_GRange)$ID) %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out"),
         up.down=ifelse(mcf7.lfc < 0,"down","up")) %>% 
  ggplot(.,aes(mcf7.padj,color=hub.io))+
  scale_color_brewer(palette="Set1")+
  geom_density()+
  facet_grid(up.down~.)
gg_tmp
ggsave("~/Documents/multires_bhicect/Poster/img/pval_dens_MCF7.svg",gg_tmp)

gg_tmp<-res_dge_tbl %>% 
  filter(ID %in% mcols(ok_dge_GRange)$ID) %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  ggplot(.,aes(mcf7.lfc,color=hub.io))+
  scale_color_brewer(palette="Set1")+
  theme_minimal()+
  geom_density()
gg_tmp
ggsave("~/Documents/multires_bhicect/weeklies/weekly61/img/lfc_dens_MCF7_2tag_edgeR.svg",gg_tmp)

#-------------------------------------
res_dge_tbl %>% 
  ggplot(.,aes(mcf7.pval))+
  geom_histogram(bins=100)+
  theme_minimal()
ggsave("~/Documents/multires_bhicect/weeklies/weekly60/img/TSS_DSEQ2_pval_hist.svg")

in_vec<-res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="in") %>% 
  dplyr::select(mcf7.lfc) %>% 
  unlist

out_vec<-res_dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="out") %>% 
  dplyr::select(mcf7.lfc) %>% 
  unlist

wilcox.test(in_vec,out_vec,alternative = "less")$p.value
