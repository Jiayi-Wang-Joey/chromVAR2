#args <- list(list.files("/mnt/plger/jwang/data/sim/01-weight", "^weight-.*", full.names=TRUE), "plts/wgt-corr.pdf")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(SummarizedExperiment)
  library(reshape2)
})

# read data
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
#se <- lapply(fs, readRDS)
info <- lapply(fs, \(x) unlist(strsplit(basename(x), "-"))) 
info <- data.frame(do.call(rbind, info))
splt <- strsplit(as.character(info$X2), ",")
info$motif <- sapply(splt, `[`, 1)
info$effect <- sapply(splt, `[`, 2)
info$X2 <- NULL
colnames(info) <- c("mode", "smooth", "peakWeight", "motif", "effect")
info$peakWeight <- gsub("\\.rds$", "", info$peakWeight)
se <- lapply(seq_len(length(fs)), \(i) {
  if (info[i,"smooth"]=="none") {
    readRDS(fs[[i]])
  } else {
    NULL
  }
})
se <- se[!vapply(se, is.null, logical(1))]
info <- info[info$smooth=="none",]

.rename <- \(y) {
  c(sapply(seq_len(ncol(y)/2), \(x) paste0("ctrl", x)),
    sapply(seq_len(ncol(y)/2), \(x) paste0("treat", x)))
}

.rmDiag <- \(x, method = "spearman") {
  colnames(x) <- .rename(x)
  corr <- cor(x, use="pairwise", method = method)
  corr[lower.tri(corr, diag=TRUE)] <- NA
  na.omit(melt(corr))
}

df <- lapply(seq_len(length(se)), \(i) {
  x <- se[[i]]
  me <- paste0(info[i,"motif"], ",", info[i,"effect"])
  total <- readRDS(paste0("/mnt/plger/jwang/data/sim/01-total/total-", 
    me, ",peaks.rds"))
  z <- assay(total, "counts")
  y <- assay(x, "counts")
  rbind(data.frame(.rmDiag(y, method="pearson"),
    cor = "Pearson", peakWeight = info[i, "peakWeight"], 
    motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(y, method="spearman"), 
      cor = "Spearman", peakWeight = info[i, "peakWeight"], 
      motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(z, method="pearson"),
      cor = "Pearson", peakWeight = "origin", 
      motif = info[i,"motif"], effect=info[i,"effect"]),
    data.frame(.rmDiag(z, method="spearman"),
      cor = "Spearman", peakWeight = "origin", 
      motif = info[i,"motif"], effect=info[i,"effect"]))
}) %>% bind_rows()

df <- df %>%
  mutate(type = case_when(
    grepl("ctrl", Var1) & grepl("ctrl", Var2) ~ "intra",
    grepl("treat", Var1) & grepl("treat", Var2) ~ "intra",
    grepl("ctrl", Var1) & grepl("treat", Var2) ~ "inter",
    grepl("treat", Var1) & grepl("ctrl", Var2) ~ "inter",
    TRUE ~ NA_character_
  ))

ps <- lapply(split(df, df$motif), \(fd) {
  lapply(split(fd, fd$effect), \(d) 
    ggplot(d, aes(x=peakWeight, y=value, col=peakWeight)) +
      geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
      geom_boxplot(width = 0.05, 
        position = position_dodge(width = 0.8), alpha = 0.1) +
      facet_grid(cor~type, scales="free") + ggtitle(d$effect[1]) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ) |> wrap_plots(n=2) + plot_layout(guides = "collect") + plot_annotation(fd$motif[1])
})


pdf(args[[2]], onefile=TRUE, width=10, height=8)
for (p in ps) print(p); dev.off()