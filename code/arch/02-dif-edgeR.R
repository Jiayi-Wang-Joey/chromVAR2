suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(chromVAR)
  library(limma)
  library(edgeR)
})

fun <- function(x, genome, motif) {
  x <- filterPeaks(x, non_overlapping = TRUE)
  motif_ix <- matchMotifs(motif, x, genome = genome)
  y <- t(assay(motif_ix))%*%assay(x, "counts")
  group_id <- substr(colnames(y), 1, nchar(colnames(y)) - 5)
  design <- model.matrix(~ group_id)
  y <- calcNormFactors(DGEList(y))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  lrt <- glmQLFTest(fit)
  res <- topTags(lrt, n=Inf)$table
  ids <- match(rownames(res), rownames(y))
  res$name <- rownames(x)[ids]
  res$rank <- seq_len(nrow(res))
  res 
  
}