library(renv)

renv::init()
renv::install("tidyverse")
renv::install("bioc::DESeq2")
renv::install("bioc::GenomicRanges")
renv::install("bioc::fgsea")
renv::snapshot()
