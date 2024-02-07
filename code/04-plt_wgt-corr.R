#args <- list(list.files("/mnt/plger/jwang/data/01-weight", "^weight-.*", full.names=TRUE), "plts/wgt-corr.pdf")

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
se <- lapply(fs, readRDS)
info <- lapply(fs, \(x) unlist(strsplit(basename(x), "-"))) 
info <- data.frame(do.call(rbind, info))
colnames(info) <- c("mode", "motif", "row", "smooth", "peakWeight")
info$peakWeight <- gsub("\\.rds$", "", info$peakWeight)


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
  m <- info[i,"motif"]
  total <- readRDS(paste0("/mnt/plger/jwang/data/dat/01-total/total-", m, "-peaks.rds"))
  z <- assay(total, "counts")
  y <- assay(x, "counts")
  rbind(data.frame(.rmDiag(y, method="pearson"),
    cor = "Pearson", peakWeight = info[i, "peakWeight"], 
    motif = m),
    data.frame(.rmDiag(y, method="spearman"), 
      cor = "Spearman", peakWeight = info[i, "peakWeight"], 
      motif = m),
    data.frame(.rmDiag(z, method="pearson"),
      cor = "Pearson", peakWeight = "origin", 
      motif = m),
    data.frame(.rmDiag(z, method="spearman"),
      cor = "Spearman", peakWeight = "origin", 
      motif = m))
}) %>% bind_rows()

df <- df %>%
  mutate(type = case_when(
    grepl("ctrl", Var1) & grepl("ctrl", Var2) ~ "intra",
    grepl("treat", Var1) & grepl("treat", Var2) ~ "intra",
    grepl("ctrl", Var1) & grepl("treat", Var2) ~ "inter",
    grepl("treat", Var1) & grepl("ctrl", Var2) ~ "inter",
    TRUE ~ NA_character_
  ))

gg <- lapply(split(df, df$motif), \(fd) 
  ggplot(fd, aes(x=peakWeight, y=value, col=peakWeight)) +
    geom_violin(position = position_dodge(width = 0.8), trim = FALSE) + 
    geom_boxplot(width = 0.05, position = position_dodge(width = 0.8), alpha = 0.1) +
    facet_wrap(cor~type, scales="free") + ggtitle(fd$motif[1]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
) |> wrap_plots(n=2) + plot_layout(guides = "collect")

ggsave(args[[2]], gg, width=50, height=40, units="cm")




