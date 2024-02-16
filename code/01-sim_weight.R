suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(R.utils)
  library(SummarizedExperiment)
  source("~/chromVAR/R/getCounts.R")
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(rtracklayer)
  library(fields)
})

# read wcs
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mtf
e <- args$eft
smt <- strsplit(args$smt, "_")[[1]]
pkw <- args$pkw
atacFrag <- readRDS(args$frg)
rowType <- args$row
if (smt[1]=="none") aRange <- 0 else aRange <- as.numeric(smt[2])
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
    mode = "weight",
    genome = BSgenome.Hsapiens.UCSC.hg38,
    species = "Homo sapiens",
    width = 300,
    smooth = smt[1],
    nGCBins = 10,
    nWidthBins = 30,
    aRange = aRange,
    peakWeight = pkw)

saveRDS(se, args$res)