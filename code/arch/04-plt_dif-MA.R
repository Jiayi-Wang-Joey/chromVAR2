#args <- list(list.files("outs", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
  df <- readRDS(x)
  exp <- ifelse(any("AveExpr" %in% colnames(df)), "AveExpr", "logCPM")
  pval <- ifelse(any("adj.P.Val" %in% colnames(df)), "adj.P.Val", "FDR")
  df <- df[,c("logFC", exp, pval, "mode", "dif", "test", "smooth")]
  names(df)[1:3] <- c("log2FoldChange", "baseMeanLog2", "padj")
  df
})

res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)
df$method <- paste0(df$mode,",", df$smooth,",", df$dif)

ps <- lapply(split(df, df$test), \(fd) {
  lapply(split(fd, fd$method), \(d) 
         # ggmaplot(d, fdr = 0.05, fc = 2, size = 0.4,
         #          palette = c("#B31B21", "#1465AC", "darkgray"),
         #          genenames = as.vector(df$name),
         #          legend = "top", top = 20,
         #          font.label = c("bold", 11),
         #          font.legend = "bold",
         #          font.main = "bold",
         #          ggtheme = ggplot2::theme_minimal()) +
        ggplot(weightLV, aes(x=baseMeanLog2, y=log2FoldChange)) +
          geom_pointdensity() +
          geom_smooth(method="lm") +
           ggtitle(d$method[[1]]) +
           theme(plot.title = element_text(size = 7))) |> wrap_plots(ncol = 4) + 
    plot_annotation(title = fd$test[1])
})

pdf(args[[2]], onefile=TRUE, width=15, height=20)
for (p in ps) print(p); dev.off()




