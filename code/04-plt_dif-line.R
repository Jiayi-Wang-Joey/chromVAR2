#args <- list(list.files("outs", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
fs <- strsplit(args[[1]], ";")[[1]]
res <- lapply(fs, function(x) {
    df <- readRDS(x)
    test <- ifelse(df$test[1]!="MYC", df$test[1], "MAX")
    #test <- df$test[1]
    df <- df[row.names(df)==test,]
    df <- df[,c("t", "rank", "mode", "dif", "test", "smooth", "peakWeight")]
    #df <- df[,c("rank", "mode", "dif", "test", "smooth")]
    
})
res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]
df <- do.call(rbind, res)

ps <- lapply(split(df, df$test), \(fd) {
    #bl <- fd[fd$mode=="total"]
    #fd <- fd[fd$mode=="weight",]
    lapply(split(fd, fd$dif), \(d) {
        f <- d[d$mode=="weight",]
        bt <- d[d$mode=="total", "t"]
        br <- d[d$mode=="total", "rank"]
        #f <- na.omit(f)
        p1 <- ggplot(f, aes(x=as.factor(peakWeight), y=t)) + 
            geom_line(stat = "identity", alpha=0.8) +
            geom_label_repel(aes(label = round(t, 2))) +
            geom_hline(yintercept=bt, linetype="dashed", color = "blue") +
            ggtitle(f$dif[1]) 
        p2 <- ggplot(f, aes(x=as.factor(peakWeight), y=rank)) + 
            geom_line(stat = "identity", alpha=0.8) +
            geom_label_repel(aes(label = rank)) +
            geom_hline(yintercept=br, linetype="dashed", color = "blue") 
        p1+p2
    }) |> wrap_plots(ncol = 1) + plot_annotation(fd$test[[1]])
    
})
pdf(args[[2]], onefile=TRUE, width=10, height=18)
for (p in ps) print(p); dev.off()


