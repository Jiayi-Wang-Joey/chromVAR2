suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(limma)
    library(chromVAR)
    library(edgeR)
    library(motifmatchr)
})

fun <- function(x, genome, motif) {
    x <- x[which(!is.infinite(rowSums(assay(x, "type_1")))),]
    x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome = genome)
    y <- t(assay(motif_ix))%*%assay(x, "type_1")
    group_id <- substr(colnames(y), 1, nchar(colnames(y)) - 5)
    design <- model.matrix(~ group_id)
    y <- calcNormFactors(DGEList(y))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf)
    ids <- match(rownames(res), rownames(y))
    res$rank <- seq_len(nrow(res))
    res 
    
}