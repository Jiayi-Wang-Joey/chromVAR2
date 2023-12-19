suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(R.utils)
    library(SummarizedExperiment)
    source("~/chromVAR/R/getCounts.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(rtracklayer)
    library(fields)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mot
atacFrag <- readRDS(args$frg)
rowType <- args$row
smt <- strsplit(args$smt, "_")[[1]]
mice <- c("BANP", "NR1H3", "NR1H4")
pf <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/",m,"/peaks"),
    full.names = TRUE)
peaks <- import.bed(con=pf)
if(smt[1]=="none") aRange <- 0 else aRange <- as.numeric(smt[2])
if (m %in% mice) {
    se <- getCounts(atacFrag = atacFrag,
        ranges = peaks,
        rowType = rowType,
        mode = "weight",
        genome = BSgenome.Mmusculus.UCSC.mm10,
        species = "Mus_musculus",
        width = 300,
        smooth = smt[1],
        aRange = aRange)
} else {
    se <- getCounts(atacFrag = atacFrag,
        ranges = peaks,
        rowType = rowType,
        mode = "weight",
        genome = BSgenome.Hsapiens.UCSC.hg38,
        species = "Homo sapiens",
        width = 300,
        smooth = smt[1],
        aRange = aRange)
}

saveRDS(se, args$res)