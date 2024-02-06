suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(R.utils)
    library(SummarizedExperiment)
    source("~/chromVAR/R/getCounts.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(rtracklayer)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mtf
e <- args$eft
atacFrag <- readRDS(args$frg)
rowType <- args$row
folder <- paste0(m, "_haploinsufficiency_", e, "_FALSE/peaks")
if (e !="0") {
  pf <- list.files(paste0("/mnt/plger/jwang/sim_data_es/",folder), 
    "^merged*",
    full.names = TRUE)
} else {
  pf <- "/mnt/plger/esonder/R/tf_activity_benchmark/DTFAB/simulations/data/peaks/merged_summits.bed"
}
peaks <- import.bed(con=pf)
se <- getCounts(atacFrag = atacFrag,
  ranges = peaks,
  rowType = rowType,
  mode = "total",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "Homo sapiens",
  width = 300)

saveRDS(se, args$res)