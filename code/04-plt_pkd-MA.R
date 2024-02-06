#args <- list(list.files("outs", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
    df <- readRDS(x)
    exp <- ifelse(any("AveExpr" %in% colnames(df)), "AveExpr", "logCPM")
    pval <- ifelse(any("adj.P.Val" %in% colnames(df)), "adj.P.Val", "FDR")
    df <- df[,c("logFC", exp, pval, "mode", "dif", "test", "smooth", "peakWeight")]
    names(df)[1:3] <- c("log2FoldChange", "baseMeanLog2", "padj")
    df
})

res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)
df$method <- paste0(df$mode,",", df$peakWeight,",", df$dif)

ps <- lapply(split(df, df$test), \(fd) {
        ggplot(fd, aes(x=baseMeanLog2, y=log2FoldChange)) +
            #geom_pointdensity() +
            geom_point_rast(alpha=0.5) +
            geom_smooth(method="lm") +
            ggtitle(fd$test[[1]]) +
            theme(plot.title = element_text(size = 7)) +
            facet_wrap(peakWeight~mode, scales="free") 
}) 
# |> wrap_plots(ncol = 4) + 
#     plot_annotation(title = fd$test[1])

pdf(args[[2]], onefile=TRUE, width=15, height=10)
for (p in ps) print(p); dev.off()