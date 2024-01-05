#args <- list(list.files("outs", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
    df <- readRDS(x)
    df <- df[row.names(df)==df$test[1],]
    df <- df[,c("t", "rank", "mode", "dif", "test", "smooth")]
   
})
res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)
df$method <- paste0(df$mode,",", df$smooth,",", df$dif)

pr <- ggplot(df, 
    aes(reorder(test,-sqrt(rank)), 
        reorder(method,-sqrt(rank)), fill=rank)) +
    geom_tile(col="white") +
    geom_text(aes(label = rank), size = 1.5) +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        n.breaks=3, direction=-1) +
    labs(x="motifs", y="method") +
    coord_fixed(expand=FALSE) + 
    theme_minimal(6) + theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank")

pt <- ggplot(df, 
    aes(reorder(test,abs(t)), 
        reorder(method,abs(t)), fill=t)) +
    geom_tile(col="white") +
    geom_text(aes(label = round(t, 2)), size = 1.5) +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        n.breaks=3, direction=-1) +
    labs(x="motifs", y="method") +
    coord_fixed(expand=FALSE) + 
    theme_minimal(6) + theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("t-value")

gg <- pt + pr + plot_annotation(tag_levels="a")
ggsave(args[[2]], gg, width=20, height=10, units="cm")

