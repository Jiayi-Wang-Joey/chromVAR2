#args <- list(list.files("outs/sim", "^dif-.*", full.names=TRUE), "plts/sim/dif-heatmap.pdf")
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

dfw <- df[df$mode=="weight",]
ps <- lapply(split(dfw, dfw$test), \(fd) {
  ggplot(fd, aes(effect, rank, col=peakWeight)) +
    geom_line(stat = "identity", alpha=0.8, aes(group=peakWeight)) +
    geom_point(stat = "identity", alpha=0.8, size = 1) +
    geom_label_repel(aes(label = rank)) +
    facet_grid(smooth~dif, scales = "free") + 
    theme_bw() +
    ggtitle(fd$test[1]) +
    scale_color_brewer(palette = "Set1") & theme(
      plot.margin=margin(),
      panel.border=element_rect(fill=NA),
      legend.key.size=unit(0.25, "lines")) 

})

pdf(args[[2]], onefile=TRUE, width=10, height=8)
for (p in ps) print(p); dev.off()