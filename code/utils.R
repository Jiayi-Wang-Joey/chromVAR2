
.chromVAR_bulk <- \(x, motifs) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x, 
        genome = BSgenome.Hsapiens.UCSC.hg38)
    bg <- getBackgroundPeaks(object = x, niterations = 1000)
    motif_ix <- matchMotifs(motifs, 
        x, 
        genome = BSgenome.Hsapiens.UCSC.hg38)
    dev <- computeDeviations(object = x, 
        annotations = motif_ix,
        background_peaks = bg
        )
}