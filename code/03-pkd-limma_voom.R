suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
})

fun <- function(x, genome, motif) {
    x <- x[which(!is.infinite(rowSums(counts(x)))),]
    counts <- counts(x)
    group_id <- rep(LETTERS[1:2],each=ncol(counts)/2)
    design <- model.matrix(~ group_id)
    rownames(counts) <- seq_len(nrow(counts))
    y <- calcNormFactors(DGEList(counts))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf, sort.by = "none")
    res 
}


