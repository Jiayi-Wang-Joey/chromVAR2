#args <- list(list.files("outs/sim", "^dif-.*", full.names=TRUE), "plts/sim/dif-bubble.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
})


args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
    df <- readRDS(x)
    test <- ifelse(df$test[1]!="ZNF143", df$test[1], "ZN143")
    df <- df[row.names(df)==test,]
    df <- df[,c("t", "rank", "mode", "dif", "effect", 
        "test", "smooth", "peakWeight", "adj.P.Val")]
})
res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)
df$method <- paste0(df$mode,",", df$smooth, ",", df$peakWeight,",", df$dif)
df$significance <- ifelse(df$adj.P.Val < 0.05, TRUE, FALSE)



ps <- lapply(split(df, df$test), \(fd) {
    rng <- range(fd$rank, na.rm=TRUE)
    rng <- c(
        floor(rng[1]*10)/10, 
        ceiling(rng[2]*10)/10)
    ggplot(fd, aes(x = effect, y = reorder(method, significance), 
        size = significance, color = rank)) +
        geom_point() +
        scale_color_gradientn(
            limits=rng, breaks=rng, na.value="lightgrey", 
            colors=c("black", "firebrick", "red", "pink", "moccasin")) +
        theme_minimal() + 
        ggtitle(fd$test[1]) +
        ylab("methods")
}) |> wrap_plots(ncol = 2)

ggsave(args[[2]], ps, width=30, height=20, units="cm")
# pdf(args[[2]], onefile=TRUE, width=6, height=4)
# for (p in ps) print(p); dev.off()