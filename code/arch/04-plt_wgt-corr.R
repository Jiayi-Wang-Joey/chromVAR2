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

df <- lapply(seq_len(length(se)), \(i) {
  x <- se[[i]]
  y <- assay(x, "counts")
  rename <- \(y) {
    c(sapply(seq_len(ncol(y)/2), \(x) paste0("ctrl", x)),
      sapply(seq_len(ncol(y)/2), \(x) paste0("treat", x)))
  }
  colnames(y) <- rename(y)
  m <- info[i,"motif"]
  total <- readRDS(paste0("/mnt/plger/jwang/data/01-total/total-", m, "-peaks.rds"))
  z <- assay(total, "counts")
  colnames(z) <- rename(z)
  rbind(data.frame(melt(lower.tri(cor(y, use="pairwise"), diag = FALSE)),
    cor = "Pearson", peakWeight = info[i, "peakWeight"], 
    motif = m),
    data.frame(melt(lower.tri(cor(y, use="pairwise", method = "spearman"),
      diag = FALSE)), 
      cor = "Spearman", peakWeight = info[i, "peakWeight"], 
      motif = m),
    data.frame(melt(lower.tri(cor(z, use="pairwise"), diag = FALSE)),
      cor = "Pearson", peakWeight = "origin", 
      motif = m),
    data.frame(melt(lower.tri(cor(z, use="pairwise", method="spearman"), 
      diag = FALSE)),
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




