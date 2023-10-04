
ann <- readRDS("../data/01_sim/rowAnnotations.rds")
addAnnotations <- \(motif_ix) {
    mc <- ann[ann$has_chip_motif == TRUE, ]
    mmMat <- assays(motif_ix)$motifMatches
    as(as(as(as.matrix(mmMat), "dMatrix"), "generalMatrix"), "TsparseMatrix")
    #mmMat <- as(as.matrix(mmMat), "dgTMatrix")
    mmMat <- as(as(as(as.matrix(mmMat), "dMatrix"), "generalMatrix"), "TsparseMatrix")
    mmDt <- data.table(i = mmMat@i,j = mmMat@j,x = mmMat@x)
    ann <- as.data.table(ann)
    ann[,peak_id := 1:nrow(ann)]
    mmAnnDt <- merge(mmDt, ann, all.x = TRUE, by.x = c("i"), by.y = c("peak_id"))
    mmChipDt <- subset(mmAnnDt, has_chip)
    mmNoChIPDt <- subset(mmAnnDt, !has_chip)
    mmChIPMat <- sparseMatrix(
        i = mmChipDt$i,j = mmChipDt$j+1, 
        x = rep(TRUE, nrow(mmChipDt)), 
        dims = dim(mmMat))
    mmNoChIPMat <- sparseMatrix(i = mmNoChIPDt$i,
        j = mmNoChIPDt$j+1, 
        x = rep(TRUE, nrow(mmNoChIPDt)), 
        dims = dim(mmMat))
    colnames(mmNoChIPMat) <- paste0("no_chip_", colnames(mmMat))
    colnames(mmChIPMat) <- paste0("with_chip_", colnames(mmMat))
    mmCombMat <- cbind(mmNoChIPMat, mmChIPMat)
    motif_ix_comb <- SummarizedExperiment(
        assays = SimpleList(motifMatches = mmCombMat),
        rowRanges = rowRanges(motif_ix),
        colData = data.frame(name = colnames(mmCombMat)))
}


.chromVAR_bulk <- \(x, motifs, genome, annotation) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x, 
        genome = genome)
    bg <- getBackgroundPeaks(object = x, niterations = 1500)
    motif_ix <- matchMotifs(motifs, 
        x, 
        genome = genome)
    if (annotation == TRUE) motif_ix <- addAnnotations(motif_ix)
    dev <- computeDeviations(object = x, 
        annotations = motif_ix,
        background_peaks = bg
        )
}

addChIP <- \(x) {
    nid <- grep("no_chip", x)
    chip <- rep("with_chip", length(x))
    chip[nid] <- "no_chip"
    return(chip)
}


getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
    species=c("Hsapiens","Mmusculus")){
    species <- match.arg(species)
    motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
    pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
    modf <- data.frame(row.names=names(motifs),
        TF=gsub(pat,"",names(motifs)),
        grade=gsub(".+\\.","",names(motifs)))
    modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
    modf <- modf[!duplicated(modf$TF),]
    motifs <- motifs[row.names(modf)]
    switch(match.arg(format),
        universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
        PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
            universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
            modf$TF)),
        PWMatrix=do.call(TFBSTools::PWMatrixList, 
            setNames(universalmotif::convert_motifs(motifs, 
                class="TFBSTools-PWMatrix"), modf$TF))
    )
}