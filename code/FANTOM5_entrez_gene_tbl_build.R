library(GenomicRanges)
library(vroom)
library(tidyverse)
library(DESeq2)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
out_file<-"./data/FANTOM5_entrez_gene_tbl.Rda"

FANTOM5_tpm_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_counts_ann.osc.txt",delim="\t",comment = "#",col_select = 1:7)

FANTOM5_entrez_tbl<-FANTOM5_tpm_tbl %>% 
  dplyr::select(`00Annotation`,entrezgene_id) %>% 
  dplyr::rename(peak=`00Annotation`) %>% 
  mutate(entrez.l=str_split(entrezgene_id,",")) %>% 
  mutate(entrez.id= map(entrez.l,function(x){
    str_extract(x,"[0-9]+")
  })) %>% 
  dplyr::select(peak,entrez.id) %>% 
  mutate(chr=str_split_fixed(peak,":|\\.\\.|,",4)[,1],
         start=as.numeric(str_split_fixed(peak,":|\\.\\.|,",4)[,2]),
         end=as.numeric(str_split_fixed(peak,":|\\.\\.|,",4)[,3])) %>% 
  distinct() %>% 
  filter(!(is.na(start)))
save(FANTOM5_entrez_tbl,file=out_file)
#-------------------------------------------------------------------
cage_Grange_fn<-function(cage_hmec_a){
  cage_a<-cage_hmec_a%>%filter(!(is.na(start)))
  full_cage_Grange<-   GRanges(seqnames=cage_a$chr,
                               ranges = IRanges(start=cage_a$start,
                                                end=cage_a$end,
                                                names=paste(cage_a$chr,1:nrow(cage_a),sep='_')
                               ))
  mcols(full_cage_Grange)<-tibble(entrez=cage_a$entrez.id)
  return(full_cage_Grange)
}
#-------------------------------------------------------------------

dge_file<-"./data/HMEC_MCF7_MDA_DSeq2.Rda"
out_file<-"./data/CAGE_DGE_entrez_gene_tbl.Rda"

cage_tbl<-get(base::load(dge_file))
tmp_obj<-names(mget(load(dge_file)))
rm(list=tmp_obj)
rm(tmp_obj)

res_mcf7 <- results( cage_tbl ,name = "cell.line_MCF7_vs_HMEC")
res_mda <- results( cage_tbl ,name = "cell.line_MDA_vs_HMEC")
res_mcf7_tbl<-as_tibble(res_mcf7)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mcf7))%>%dplyr::rename(mcf7.lfc=log2FoldChange,mcf7.lfc.se=lfcSE,mcf7.pval=pvalue,mcf7.padj=padj)
res_mda_tbl<-as_tibble(res_mda)%>%dplyr::select(log2FoldChange,lfcSE, pvalue, padj)%>%mutate(ID=rownames(res_mda))%>%dplyr::rename(mda.lfc=log2FoldChange,mda.lfc.se=lfcSE,mda.pval=pvalue,mda.padj=padj)
res_dge_tbl<-res_mcf7_tbl%>%full_join(.,res_mda_tbl)
rm(res_mcf7_tbl,res_mda_tbl,res_mda,res_mcf7)
rm(cage_tbl)

entrez_dge_tbl<-res_dge_tbl %>% 
  left_join(.,FANTOM5_entrez_tbl,by=c("ID"="peak")) %>% 
  dplyr::select(ID,mcf7.lfc,mcf7.padj,mda.lfc,mda.padj,entrez.id) %>% 
  mutate(chr=str_split_fixed(ID,":|\\.\\.|,",4)[,1],
         start=as.numeric(str_split_fixed(ID,":|\\.\\.|,",4)[,2]),
         end=as.numeric(str_split_fixed(ID,":|\\.\\.|,",4)[,3]))


save(entrez_dge_tbl,file=out_file)
