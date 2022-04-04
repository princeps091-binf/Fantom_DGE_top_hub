# # cell-line DGE analysis
library(DESeq2)
library(tidyverse)
library(viridis)
library(kableExtra)
library(fgsea)
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
save(dds,file = "~/Documents/multires_bhicect/DGE_HiC_integration/HMEC_MCF7_MDA_DSeq2.Rda")
resultsNames(dds)
res_mcf7 <- results( dds ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( dds ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(col_data_tbl,count_data_tbl,ddsFullCountTable)

res_dge_tbl %>% arrange(desc(mcf7.lfc)) %>% mutate(idx=1:n(),signif=ifelse(is.na(mcf7.padj),"out",ifelse(mcf7.padj<0.01,"FDR < 0.01","FDR > 0.01"))) %>% 
  filter(sign(mcf7.lfc+mcf7.lfc.se)*sign(mcf7.lfc-mcf7.lfc.se)==1) %>% 
  ggplot(aes(x=mcf7.lfc-mcf7.lfc.se,xend=mcf7.lfc+mcf7.lfc.se,y=idx,yend=idx,color=signif))+
  geom_segment(alpha=1,size=0.05)+
  geom_vline(xintercept = 0,linetype=2)+facet_wrap(signif~.,scales="free")+
  theme(legend.position = "none")

res_dge_tbl %>% arrange(desc(mda.lfc+mda.lfc.se)) %>% 
  mutate(idx=1:n(),signif=ifelse(is.na(mda.padj),"out",ifelse(mda.padj<0.01,"FDR < 0.01","FDR > 0.01"))) %>% 
  #filter(!(is.na(mcf7.padj))) %>% 
  ggplot(aes(x=mda.lfc-mda.lfc.se,xend=mda.lfc+mda.lfc.se,y=idx,yend=idx,color=signif))+
  geom_segment(alpha=0.5,size=0.05)+
  geom_vline(xintercept = 0,linetype=2)+facet_wrap(signif~.)+
  theme(legend.position = "none")

#------------
# Logistic regression to evaluate the extent of LFC agreement as a function of DGE significance
logit_dat<-res_dge_tbl%>%mutate(agg=sign(mcf7.lfc)*sign(mda.lfc))%>%mutate(agg=ifelse(agg<0,0,1))%>%dplyr::select(ID,agg,mcf7.pval,mcf7.padj,mda.pval,mda.padj)%>%filter(!(is.na(mcf7.padj))&!(is.na(mda.padj)))

mylogit <- glm(agg ~ mcf7.padj + mda.padj , data = logit_dat, family = "binomial")

summary(mylogit)
confint(mylogit)
exp(coef(mylogit))
expand_grid(mcf7.padj=seq(1e-3,0.1,length.out=100),mda.padj=seq(1e-2,0.1,length.out=100))
newdata1 <- expand_grid(mcf7.padj=seq(1e-2,0.1,length.out=100),mda.padj=seq(1e-2,0.1,length.out=100))
newdata1<-newdata1%>%mutate(agg.pred=predict(mylogit, newdata = newdata1, type = "response"))
newdata1%>%ggplot(.,aes(mcf7.padj,mda.padj,fill=agg.pred))+geom_raster()

gg_lfc<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc))+geom_point(alpha=0.1)
ggsave("~/Documents/multires_bhicect/weeklies/weekly37/img/mcf_mda_lfc_point.png",gg_lfc)
gg_cor<-do.call(bind_rows,lapply(seq(1e-10,1,length.out=500),function(x){
  tmp_tbl<-res_dge_tbl%>%filter(mcf7.padj<=x | mda.padj<= x)
  tibble(cor=cor(tmp_tbl$mcf7.lfc,tmp_tbl$mda.lfc),padj=x)
}))%>%ggplot(.,aes(padj,cor))+geom_line()
ggsave("~/Documents/multires_bhicect/weeklies/weekly37/img/mcf_mda_cor.png",gg_cor)

gg_bar<-do.call(bind_rows,lapply(seq(1e-10,1,length.out=500),function(x){
  tmp_tbl<-res_dge_tbl%>%filter(mcf7.padj<=x | mda.padj<= x)
  return(tmp_tbl%>%mutate(agg=sign(mcf7.lfc)*sign(mda.lfc))%>%group_by(agg)%>%summarise(n=n())%>%mutate(padj=x))
}))%>%ggplot(.,aes(padj,n,fill=as.factor(agg)))+geom_bar(stat="identity",position='fill')
ggsave("~/Documents/multires_bhicect/weeklies/weekly37/img/mcf_mda_bar.png",gg_bar)

# In vs Out CAGE-enriched clusters

load('~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/dagger_mres_fdr01_multi_cage_multi_bin_peak_shuffle_topcl_tbl.Rda')
res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
#Compute CAGE peak composition of each top cluster
res_dge_tbl<-res_dge_tbl%>%mutate(chr=unlist(lapply(strsplit(.$ID,split=":"),'[',1)))
tmp_coord<-unlist(lapply(strsplit(unlist(lapply(strsplit(res_dge_tbl$ID,split=":"),'[',2)),split=','),'[',1))
tmp_start<-as.numeric(unlist(lapply(strsplit(tmp_coord,split="\\.."),'[',1)))
tmp_end<-as.numeric(unlist(lapply(strsplit(tmp_coord,split="\\.."),'[',2)))
res_dge_tbl<-res_dge_tbl%>%mutate(start=tmp_start,end=tmp_end)
rm(tmp_coord,tmp_start,tmp_end)

cage_set_tbl<-res_dge_tbl%>%filter(chr %in% unique(top_cl_tbl$chr))
cage_dge_Grange<-   GRanges(seqnames=cage_set_tbl$chr,
                            ranges = IRanges(start=cage_set_tbl$start,
                                             end=cage_set_tbl$end
                            ))
#mcols(cage_dge_Grange)<-tibble(lfc=cage_set_tbl$log2FoldChange,p.score=-log10(cage_set_tbl$pvalue),mag=cage_set_tbl$baseMean)%>%mutate(dge.prank=percent_rank(p.score),mag.prank=percent_rank(mag))
mcols(cage_dge_Grange)<-cage_set_tbl%>%dplyr::select(mcf7.lfc, mcf7.pval, mcf7.padj, ID, mda.lfc, mda.pval, mda.padj)

## Produce for each top cluster a corresponding GRange object
chr_cl_l<-GRangesList()
for(chromo in unique(top_cl_tbl$chr)){
  print(chromo)
  load(paste0(res_file,chromo,"_spec_res.Rda"))
  tmp_cl<- unlist(top_cl_tbl%>%filter(chr==chromo)%>%dplyr::select(node)) 
  
  tmp_res<-unlist(lapply(strsplit(tmp_cl,split="_"),'[',1))
  chr_cl_tbl<-tibble(chr=chromo,res=tmp_res,cl=tmp_cl,bins=lapply(chr_spec_res$cl_member[tmp_cl],as.numeric))
  
  
  #Convert to Grange objec
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl,c("chr_cl_tbl","res_num"))
  res_list<-GRangesList(parLapply(cl,1:nrow(chr_cl_tbl),function(x){
    
    cl_Grange<-   GRanges(seqnames=chr_cl_tbl$chr[x],
                          ranges = IRanges(start=chr_cl_tbl$bins[[x]],
                                           end=chr_cl_tbl$bins[[x]] + res_num[chr_cl_tbl$res[x]]-1
                          ))
    
  }))
  stopCluster(cl)
  rm(cl)
  
  chr_cl_l<-append(chr_cl_l,res_list)
  
}
rm(res_list,tmp_cl,tmp_res,chr_cl_tbl,chr_spec_res)

tmp_inter<-findOverlaps(chr_cl_l,cage_dge_Grange)
tmp_inter_tbl<-as_tibble(tmp_inter)%>%group_by(queryHits)%>%
  summarise(ID=list(cage_dge_Grange@elementMetadata$ID[subjectHits]), mcf7.lfc= list(cage_dge_Grange@elementMetadata$mcf7.lfc[subjectHits]), mda.lfc= list(cage_dge_Grange@elementMetadata$mda.lfc[subjectHits]),mcf7.padj=list(cage_dge_Grange@elementMetadata$mcf7.padj[subjectHits]),mda.padj=list(cage_dge_Grange@elementMetadata$mda.padj[subjectHits]))%>%
  mutate(chr=top_cl_tbl$chr[queryHits],cl=top_cl_tbl$node[queryHits])%>%dplyr::select(-queryHits)
tmp_inter_tbl<-tmp_inter_tbl%>%mutate(res=unlist(lapply(strsplit(.$cl,split='_'),'[',1)))

chr_cl_tbl<-tmp_inter_tbl%>%mutate(cage_n=countOverlaps(chr_cl_l,cage_dge_Grange))

chr_cl_tbl<-chr_cl_tbl%>%mutate(cage.tbl=lapply(1:nrow(chr_cl_tbl),function(x){
  return(tibble(ID=chr_cl_tbl$ID[[x]],mcf7.lfc=chr_cl_tbl$mcf7.lfc[[x]],mcf7.padj=chr_cl_tbl$mcf7.padj[[x]],mda.lfc=chr_cl_tbl$mda.lfc[[x]],mda.padj=chr_cl_tbl$mda.padj[[x]]))
}))

chr_cl_tbl<-chr_cl_tbl%>%dplyr::select(chr,res,cl,cage_n,cage.tbl)
#---------------------------------------------------------------------
# Examine enrichment of DGE signal in CAGE-enriched clusters
res_cage_ID<-chr_cl_tbl%>%unnest(cols=c(cage.tbl))%>%group_by(res)%>%distinct(ID)%>%summarise(ID=list(ID))
## Set-up GSEA analysis
## Examine if the genes/peaks contained in clusters are enriched for downregulated genes
tmp_path_l<-res_cage_ID$ID
names(tmp_path_l)<-res_cage_ID$res


## MCF7
### padj
cage_peak_tbl<-res_dge_tbl%>%filter( !(is.na(mcf7.padj)))%>%arrange(desc(-log10(mcf7.padj)))
lfc_rank<--log10(cage_peak_tbl$mcf7.padj)
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
    save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_pval_mcf7.png"))

### LFC
cage_peak_tbl<-res_dge_tbl%>%filter( !(is.na(mcf7.lfc)))%>%arrange(mcf7.lfc)
lfc_rank<-cage_peak_tbl$mcf7.lfc
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
  save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_lfc_mcf7.png"))


## MDA
### padj
cage_peak_tbl<-res_dge_tbl%>%filter( !(is.na(mda.padj)))%>%arrange(desc(-log10(mda.padj)))
lfc_rank<--log10(cage_peak_tbl$mda.padj)
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
  save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_pval_mda.png"))

### LFC
cage_peak_tbl<-res_dge_tbl%>%filter( !(is.na(mda.lfc)))%>%arrange(mda.lfc)
lfc_rank<-cage_peak_tbl$mda.lfc
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
  save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_lfc_mda.png"))
## Reminder of resulter obtained when focusing on down-regulated genes
cage_peak_tbl<-res_dge_tbl%>%filter(mcf7.lfc>0 & !(is.na(mcf7.padj)))%>%arrange(desc(-log10(mcf7.padj)))
lfc_rank<--log10(cage_peak_tbl$mcf7.padj)
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
  save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_up_pval_mcf7.png"))

cage_peak_tbl<-res_dge_tbl%>%filter(mda.lfc>0 & !(is.na(mda.padj)))%>%arrange(desc(-log10(mda.padj)))
lfc_rank<--log10(cage_peak_tbl$mda.padj)
names(lfc_rank)<-as.character(cage_peak_tbl$ID)

fgseaRes <- fgseaSimple(
  pathways=tmp_path_l,
  stats=lfc_rank,
  nperm=1e4,
  minSize = 5,
  maxSize = Inf,
  scoreType = "std",
  nproc = 3,
  gseaParam = 1,
  BPPARAM = NULL
)
tmp_tbl<-fgseaRes%>%dplyr::select(pathway,padj,NES)
kbl(tmp_tbl)%>%column_spec(3, width = "5em")%>%
  save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly38/img/in_cluster_GSEA_up_pval_mda.png"))

#---------------------------------------------------------------------
res_cage_ID<-chr_cl_tbl%>%unnest(cols=c(cage.tbl))%>%group_by(res)%>%distinct(ID)%>%summarise(ID=list(ID))
res_cage_ID_l<-res_cage_ID$ID
names(res_cage_ID_l)<-res_cage_ID$res
res_dge_tbl<-res_dge_tbl%>%
                mutate(io.1Mb=ifelse(ID %in% res_cage_ID_l[['1Mb']],"in","out"))%>%
                mutate(io.500kb=ifelse(ID %in% res_cage_ID_l[['500kb']],"in","out"))%>%
                mutate(io.100kb=ifelse(ID %in% res_cage_ID_l[['100kb']],"in","out"))%>%
                mutate(io.50kb=ifelse(ID %in% res_cage_ID_l[['50kb']],"in","out"))%>%
                mutate(io.10kb=ifelse(ID %in% res_cage_ID_l[['10kb']],"in","out"))%>%
                mutate(io.5kb=ifelse(ID %in% res_cage_ID_l[['5kb']],"in","out"))
library(cowplot)
gg_5<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.5kb))+geom_point()+facet_wrap(io.5kb~.)
gg_10<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.10kb))+geom_point()+facet_wrap(io.10kb~.)
gg_50<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.50kb))+geom_point()+facet_wrap(io.50kb~.)
gg_100<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.100kb))+geom_point()+facet_wrap(io.100kb~.)
gg_500<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.500kb))+geom_point()+facet_wrap(io.500kb~.)
gg_1<-res_dge_tbl%>%ggplot(.,aes(mcf7.lfc,mda.lfc,color=io.1Mb))+geom_point()+facet_wrap(io.1Mb~.)
plot_grid(gg_5,gg_10,gg_50,gg_100,gg_500,gg_1,nrow=3,ncol=2,labels=c("5kb","10kb","50kb","100kb","500kb","1Mb"))

res_dge_tbl%>%select(-c(mcf7.pval,mda.pval,mcf7.padj,mda.padj,chr,start,end))%>%
                pivot_longer(cols = starts_with("io"),names_to = "res",values_to = "io")%>%
                  group_by(res,io)%>%summarise(c=cor(mcf7.lfc,mda.lfc))%>%
                    ggplot(.,aes(res,c,fill=io))+geom_bar(stat="identity",position="dodge")
mres_dge_tbl<-res_dge_tbl%>%select(-c(mcf7.pval,mda.pval,chr,start,end))%>%filter(!(is.na(mcf7.padj)) & !(is.na(mda.padj)))%>%
  pivot_longer(cols = starts_with("io"),names_to = "res",values_to = "io")
gg_cor<-do.call(bind_rows,lapply(seq(1e-10,1,length.out=500),function(x){
  tmp_tbl<-mres_dge_tbl%>%filter(mcf7.padj<=x | mda.padj<= x)
  return(tmp_tbl%>%group_by(res,io)%>%summarise(cor=cor(mcf7.lfc,mda.lfc))%>%mutate(padj=x))
}))%>%mutate(res=factor(res,levels=paste0("io.",res_set)))%>%ggplot(.,aes(padj,cor,color=io))+geom_line()+facet_wrap(res~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly37/img/mcf_mda_cor_io.png",gg_cor)
#-----------------------------------------------------------------------------------------
load("~/Documents/multires_bhicect/data/epi_data/CAGE/CAGE_tss_ann_smpl.Rda")
res_dge_tbl<-res_dge_tbl%>%left_join(.,tss_ann,by=c('ID'='namess'))
rm(dds)

cage_active_genes<-cage_count_hmec%>%dplyr::select(ID)%>%
                    full_join(.,cage_count_mcf7%>%dplyr::select(ID))%>%
                      full_join(.,cage_count_MDA%>%dplyr::select(ID))%>%
                        left_join(.,tss_ann%>%dplyr::select(namess,ENSG),by=c("ID"="namess"))%>%
                          unnest(cols = c(ENSG))

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")

genes_hm <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=unique(cage_active_genes$ENSG),
  mart=ensembl)
cage_active_genes<-cage_active_genes%>%left_join(.,genes_hm,by=c("ENSG"="ensembl_gene_id"))
# Input and format gene sets of interest
hm_gene_set<-as_tibble(read.table("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/h.all.v7.3.entrez.gmt",header = F,sep = "\t",fill=T))

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
})
clusterExport(cl,c('hm_gene_set'))
hallmark_set<-parLapply(cl,1:nrow(hm_gene_set),function(x){
  tmp<-hm_gene_set[x,-c(1,2)]
  return(tmp[!(is.na(tmp))])
})
stopCluster(cl)
rm(cl)
names(hallmark_set)<-hm_gene_set$V1
hallmark_set<-lapply(hallmark_set,as.character)

# Subset the Biological Processes gene sets
GOBP_set<-hallmark_set[grep("GOBP",names(hallmark_set),value=T)]
GOBP_set<-lapply(GOBP_set,as.character)

gene_lfc<-res_dge_tbl%>%filter(!(is.na(mda.padj)))%>%
            dplyr::select(mcf7.lfc,mcf7.padj,ID,ENSG)%>%
                unnest(col=c(ENSG))%>%group_by(ENSG)%>%slice_min(mcf7.padj)%>%
                  filter(mcf7.lfc > 0 & mcf7.padj<0.01)%>%left_join(.,genes_hm,by=c("ENSG"="ensembl_gene_id"))

cl_set_gene<-as.character(unlist(gene_lfc%>%ungroup()%>%distinct(entrezgene_id)))
cage_active_genes_vec<-as.character(unlist(cage_active_genes%>%distinct(entrezgene_id)))

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
  print('node ready')
})
clusterExport(cl,c('cl_set_gene','cage_active_genes_vec'))
go_pval<-parLapply(cl,hallmark_set,function(tmp_set){
  hitInSample<-sum(cl_set_gene %in% tmp_set)
  sampleSize<-length(cl_set_gene)
  hitInPop<-sum(cage_active_genes_vec %in% tmp_set)
  failInPop<-length(cage_active_genes_vec) - hitInPop
  return(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
})
stopCluster(cl)
rm(cl)
top_set<-sort(p.adjust(unlist(go_pval),method='fdr'))
path_tbl<-tibble(Gene.Set=names(top_set),p.adjust=top_set)%>%filter(p.adjust<=0.05)
kbl(path_tbl%>%arrange(p.adjust)%>%slice_head(n=40))%>%save_kable(paste0("~/Documents/multires_bhicect/weeklies/weekly39/img/mcf7_up_FDR01_hm.png"))

