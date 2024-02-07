suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    source("code/utils.R")
    library(MotifDb)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
se <- readRDS(args$dat)
source(args$fun)

Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
res <- fun(se, 
    genome = BSgenome.Hsapiens.UCSC.hg38, 
    motif = Hmotifs)

if (is.null(args$smt)) {
    df <- data.frame(res, test=args$mtf,  
      effect=args$eft, dif=args$dif, 
      mode="total", row=args$row, 
      smooth="none", peakWeight="none",
        row.names=rownames(res))
} else {
    df <- data.frame(res, test=args$mtf, 
      effect=args$eft, dif=args$dif, 
      mode="weight", row=args$row, 
      smooth=args$smt, peakWeight=args$pkw,
      row.names=rownames(res))
}

saveRDS(df, args$res)