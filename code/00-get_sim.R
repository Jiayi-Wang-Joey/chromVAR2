suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(R.utils)
  source("~/chromVAR/R/getCounts.R")
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mtf
e <- args$eft

fd <- paste0("/mnt/plger/jwang/sim_data_es/",
    m, "_haploinsufficiency_", e, "_FALSE/seq_files")
fs <- list.files(fd, "\\.bam$|\\.bed$", full.names=TRUE)
lst <- as.list(fs)
names(lst) <- basename(fs)
x <- .importFragments(lst)
saveRDS(x, args$res)
