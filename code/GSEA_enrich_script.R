library(tidyverse)
library(GenomicRanges)
library(furrr)
library(parallel)
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
dge_tbl_file<-"./data/tss_2tag_edgeR_CAGE_DGE_entrez_gene_tbl.Rda"
#-------------------------------------------------------------------------------------------------------
dge_tbl<-tbl_in_fn(dge_tbl_file)

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


dge_Grange<-   GRanges(seqnames=dge_tbl$chr,
                       ranges = IRanges(start=dge_tbl$start,
                                        end=dge_tbl$end
                       ))
mcols(dge_Grange)<-dge_tbl %>% dplyr::select(-c(chr,start,end))

ok_dge_GRange<-dge_Grange[unique(c(which(!(is.na(mcols(dge_Grange)$mcf7.padj))),which(!(is.na(mcols(dge_Grange)$mda.padj)))))]

in_cl_peak<-mcols(ok_dge_GRange)$ID[unique(queryHits(findOverlaps(ok_dge_GRange,cl_GRange)))]


all_peak_entrez_vec<-dge_tbl %>% 
  dplyr::select(ID,entrez.id) %>% 
  unnest(cols=c(entrez.id)) %>% 
  distinct(entrez.id)%>% unlist

dge_entrez_vec<-dge_tbl %>% 
  filter(ID %in% mcols(ok_dge_GRange)$ID) %>%
  dplyr::select(ID,entrez.id) %>% 
  unnest(cols=c(entrez.id)) %>% 
  distinct(entrez.id)%>% unlist

hub_entrez_vec<-dge_tbl %>% 
  filter(ID %in% in_cl_peak) %>%
  dplyr::select(ID,entrez.id) %>% 
  unnest(cols=c(entrez.id)) %>% 
  distinct(entrez.id) %>% unlist

hub_up_entrez_vec<-dge_tbl %>% 
  filter(ID %in% in_cl_peak) %>%
  filter(!(is.na(mcf7.padj))) %>% 
  filter(mcf7.lfc > 0 ) %>% 
#  filter(!(is.na(mda.padj))) %>% 
#  filter(mda.lfc > 0 ) %>% 
  dplyr::select(ID,entrez.id) %>% 
  unnest(cols=c(entrez.id)) %>% 
  distinct(entrez.id) %>% unlist


#-------------------------------------------------------------------
GO_set_enrich_fn<-function(cl_set_gene,cage_active_genes_vec,GOBP_set){
  fn_env<-environment()
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print('node ready')
  })
  clusterExport(cl,c('cl_set_gene','cage_active_genes_vec'),envir = fn_env)
  go_pval<-parLapply(cl,GOBP_set,function(tmp_set){
    hitInSample<-sum(cl_set_gene %in% tmp_set)
    sampleSize<-length(cl_set_gene)
    hitInPop<-sum(cage_active_genes_vec %in% tmp_set)
    failInPop<-length(cage_active_genes_vec) - hitInPop
    p_val<-phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    OR_GO<-(hitInSample/sampleSize)/(hitInPop/length(cage_active_genes_vec))
    return(tibble(p.val=p_val,OR=OR_GO,in.gene=hitInSample))
  })
  stopCluster(cl)
  rm(cl)
  path_tbl<-do.call(bind_rows,go_pval)%>%mutate(Gene.Set=names(go_pval),FDR=p.adjust(p.val,method='fdr'))%>%dplyr::select(Gene.Set,FDR,OR,in.gene)
  return(path_tbl)
}

#------------------------------

gene_set_file<-"~/Documents/multires_bhicect/GO_enrichment_viz/data/Hallmark_gene_set_l.Rda"




Gene_set_l<-tbl_in_fn(gene_set_file)

full_bg_vec<-unique(unlist(Gene_set_l))

foreground_gene_vec<-hub_entrez_vec

background_gene_vec<-dge_entrez_vec


path_tbl<-GO_set_enrich_fn(foreground_gene_vec,background_gene_vec,Gene_set_l)
print(path_tbl %>% 
        filter(FDR<=0.01) %>% 
#        arrange(desc(OR))
              arrange(FDR)
      ,n=100)
#------------------------------
library(fgsea)
all_dge_peak<-mcols(ok_dge_GRange)$ID
in_cl_peak
rank_tbl<-dge_tbl %>% 
  filter(ID %in% in_cl_peak) %>% 
  dplyr::select(ID,mcf7.lfc,mcf7.padj,entrez.id) %>% 
  mutate(peak.score=sign(mcf7.lfc)*-log10(mcf7.padj)) %>% 
  unnest(cols=c(entrez.id)) %>% 
  filter(!(is.na(entrez.id))) %>% 
  group_by(entrez.id) %>% 
  summarise(entrez.score=mean(peak.score,na.rm=T)) %>% 
  filter(!(is.nan(entrez.score)))

rank_tbl<-dge_tbl %>% 
#  filter(ID %in% in_cl_peak) %>% 
  dplyr::select(ID,mda.lfc,mda.padj,entrez.id) %>% 
  mutate(peak.score=sign(mda.lfc)*-log10(mda.padj)) %>% 
  unnest(cols=c(entrez.id)) %>% 
  filter(!(is.na(entrez.id))) %>% 
  group_by(entrez.id) %>% 
  summarise(entrez.score=mean(peak.score,na.rm=T)) %>% 
  filter(!(is.nan(entrez.score)))


entrez_rank<-rank_tbl$entrez.score
names(entrez_rank)<-rank_tbl$entrez.id
fgseaRes <- fgseaSimple(
  pathways=Gene_set_l,
  stats=entrez_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
print(as_tibble(fgseaRes) %>% filter(padj<=0.01) %>% 
  arrange(NES),n=100)

topPathways <- as_tibble(fgseaRes) %>% filter(padj<=0.01) %>% 
  arrange(NES) %>% dplyr::select(pathway) %>% unlist
plotGseaTable(Gene_set_l[topPathways], entrez_rank, fgseaRes, 
              gseaParam=1)

library(formattable)
formattable(print(as_tibble(fgseaRes) %>% filter(padj<=0.01) %>% 
                    dplyr::select(pathway,padj,NES) %>% 
                    arrange(NES),n=100))


data(examplePathways)
data(exampleRanks)
fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=1000,
                  minSize=15, maxSize=100)
topPathways <- fgseaRes[head(order(pval), n=15)][order(NES), pathway]

plotGseaTable(examplePathways[topPathways], exampleRanks,
              fgseaRes, gseaParam=0.5)
