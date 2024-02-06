suppressPackageStartupMessages({
    source("code/utils.R")
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- args$mot
se <- readRDS(args$dat)
source(args$fun)

res <- fun(se)

if (is.null(args$smt)) {
    df <- data.frame(res, test=m, 
        dif=args$pkd, mode="total", 
        row=args$row, smooth="none", peakWeight="none",
        row.names=rownames(res))
} else {
    df <- data.frame(res, test=m, 
        dif=args$pkd, mode="weight", 
        row=args$row, smooth=args$smt, peakWeight=args$pkw,
        row.names=rownames(res))
}

saveRDS(df, args$res)
