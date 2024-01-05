suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(R.utils)
  source("~/chromVAR/R/getCounts.R")
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mot

fs <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/",
                          m, "/seq_files"), pattern = "\\.bam$|\\.bed$",
                   full.names = TRUE)
lst <- as.list(fs)
names(lst) <- basename(fs)
x <- .importFragments(lst)
  
saveRDS(x, args$res)
