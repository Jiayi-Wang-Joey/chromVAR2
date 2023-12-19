suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("code/utils.R")
    library(MotifDb)
})

mice <- c("BANP", "NR1H3", "NR1H4")
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mot
se <- readRDS(args$dat)
source(args$fun)

if (m %in% mice) {
    Mmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Mmusculus")
    res <- fun(se, 
        genome = BSgenome.Mmusculus.UCSC.mm10, 
        motif = Mmotifs)
} else {
    Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    res <- fun(se, 
        genome = BSgenome.Hsapiens.UCSC.hg38, 
        motif = Hmotifs)
}

if (is.null(args$smt)) {
    df <- data.frame(res, test=m, 
        dif=args$dif, mode="total", 
        row=args$row, smooth="none", 
        row.names=rownames(res))
} else {
    df <- data.frame(res, test=m, 
        dif=args$dif, mode="weight", 
        row=args$row, smooth=args$smt, 
        row.names=rownames(res))
}

saveRDS(df, args$res)
