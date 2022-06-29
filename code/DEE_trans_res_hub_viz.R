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
dge_file<-"./data/enh_HMEC_MCF7_MDA_DSeq2.Rda"
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

dge_coord_tbl<-dge_tbl %>%
  dplyr::select(ID) %>% 
  mutate(chr=str_split_fixed(ID,":|\\.\\.|,|-",4)[,1],
         start=as.numeric(str_split_fixed(ID,":|\\.\\.|,|-",4)[,2]),
         end=as.numeric(str_split_fixed(ID,":|\\.\\.|,|-",4)[,3]))
dge_Grange<-   GRanges(seqnames=dge_coord_tbl$chr,
                       ranges = IRanges(start=dge_coord_tbl$start,
                                        end=dge_coord_tbl$end
                       ))
mcols(dge_Grange)<-dge_tbl

ok_dge_GRange<-dge_Grange[unique(c(which(!(is.na(mcols(dge_Grange)$mcf7.padj))),which(!(is.na(mcols(dge_Grange)$mda.padj)))))]

in_cl_peak<-mcols(ok_dge_GRange)$ID[unique(queryHits(findOverlaps(ok_dge_GRange,cl_GRange)))]

gg_tmp<-dge_tbl %>% 
  filter(ID %in% mcols(ok_dge_GRange)$ID) %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  ggplot(.,aes(mcf7.lfc,color=hub.io))+
  scale_color_brewer(palette="Set1")+
  theme_minimal()+
  geom_density()
gg_tmp

in_vec<-dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="in") %>% 
  dplyr::select(mcf7.lfc) %>% 
  unlist

out_vec<-dge_tbl %>% 
  mutate(hub.io=ifelse(ID %in% in_cl_peak,"in","out")) %>% 
  filter(hub.io=="out") %>% 
  dplyr::select(mcf7.lfc) %>% 
  unlist

wilcox.test(in_vec,out_vec,alternative = "less")
