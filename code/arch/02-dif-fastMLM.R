# functions from Pierre-Luc
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(BiocParallel)
    library(limma)
})
fastMLM <- function(accmat, annotation, useIntercept=TRUE, poisson=FALSE){
    stopifnot(is.matrix(accmat) && is.matrix(annotation))
    if(useIntercept) annotation <- cbind(rep(1L,nrow(annotation)), annotation)
    res <- lapply(seq_len(ncol(accmat)), function(i){
        if(!isTRUE(poisson)){
            mod <- RcppArmadillo::fastLmPure(annotation, accmat[,i])
            tvals <- mod$coefficients/mod$stderr
        }else{
            if(FALSE && require("Rfast", quietly=TRUE, include.only="glm_poisson")){
                mod <- glm_poisson(a, y[,1], full=TRUE)$info
            }else{
                mod <- glm(accmat[,i]~0+annotation, family="poisson")
                mod <- coef(summary(mod))
            }
            tvals <- mod[,1]/mod[,2]
        }
        tvals
    })
    res <- matrix(unlist(res), nrow=ncol(annotation))
    row.names(res) <- colnames(annotation)
    if(useIntercept) res <- res[-1,,drop=FALSE]
    colnames(res) <- colnames(accmat)
    res
}

fun <- function(x, 
    genome, motif) {
    counts <- assay(x, "counts")
    design <- rep(LETTERS[1:2],each=ncol(counts)/2)
    #counts$condition <- design
    motif_ix <- matchMotifs(motif, x, genome = genome)
    #y <- t(assay(motif_ix))%*%assay(x, "counts")
    e <- fastMLM(counts, 
        as.matrix(assay(motif_ix)))
    fit <- eBayes(lmFit(e, model.matrix(~design)))
    res <- topTable(fit, number=Inf)
    #res <- data.frame(row.names=row.names(res), logFC=res$logFC, t=res$t,
    #    p=res$P.Value, padj=res$adj.P.Val)
    #res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    #runtime <- proc.time()-ptm
    #raw <- list(runtime=runtime, runtime2=runtime, obj1=e)
    return(res)
}
