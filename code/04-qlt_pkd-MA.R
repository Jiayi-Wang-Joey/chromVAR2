#args <- list(list.files("outs/sim", "^pkd-.*", full.names=TRUE), "plts/pkd-MA.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    library(patchwork)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
    df <- readRDS(x)
    exp <- ifelse(any("AveExpr" %in% colnames(df)), "AveExpr", "logCPM")
    pval <- ifelse(any("adj.P.Val" %in% colnames(df)), "adj.P.Val", "FDR")
    df <- df[,c("logFC", exp, pval, "mode", 
        "test", "smooth", "peakWeight", "effect")]
    names(df)[1:3] <- c("log2FoldChange", "baseMeanLog2", "padj")
    df
})

res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)

gg <- list(
    geom_point_rast(alpha=0.5),
    geom_smooth(method="lm"),
    theme(plot.title = element_text(size = 7))
)

df1 <- df[df$mode=="weight", ]
df2 <- df[df$mode=="total",]
ps <- lapply(split(df1, df1$test), \(fd) {
    dd <- df2[df2$test==fd$test[1],]
    p1 <- ggplot(fd, aes(x=baseMeanLog2, y=log2FoldChange)) +
        gg +
        facet_grid(effect~peakWeight+smooth, scales="free")
    
    p2 <- ggplot(dd, aes(x=baseMeanLog2, y=log2FoldChange)) + 
        gg + ggtitle("original counts")
    list(p1,p2) |> wrap_plots(ncol=1) + plot_annotation(fd$test[1]) +
        plot_layout(heights=c(2,1))
})


pdf(args[[2]], onefile=TRUE, width=8, height=13)
for (p in ps) print(p); dev.off()