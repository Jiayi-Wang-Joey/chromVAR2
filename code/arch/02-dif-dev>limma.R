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
    y <- cpm(calcNormFactors(DGEList(y)))
    devs <- (y-rowMeans(y))/rowMeans(y)
    fit <- eBayes(lmFit(devs, design))
    res <- topTable(fit, n = Inf)
    #ids <- match(rownames(res), rownames(y))
    #res$name <- rownames(x)$name[ids]
    res$rank <- seq_len(nrow(res))
    res 
    
}