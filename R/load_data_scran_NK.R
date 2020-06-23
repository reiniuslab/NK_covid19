## LOAD DATA
library(data.table)
library(Seurat)
fls <- list.files(pattern = ".h5$")
balf.data <- lapply(fls, Read10X_h5)
names(balf.data) <- sapply(strsplit(fls,"_"),"[[",2)

meta.cell <- fread("all.cell.annotation.meta.txt")[!grepl("^GSM",sample)]
meta.nk <- meta.cell[celltype == "NK"]
idx.nk <- split(meta.nk[,gsub("_.*","",ID)],meta.nk$sample)

nk.data <- mapply(balf.data[names(idx.nk)], idx.nk, FUN = function(x,y) x[1:min(sapply(balf.data,nrow)),gsub("-.*","",colnames(x)) %in% y, drop=F] )

library(magrittr)
library(scater)
library(scran)
sce.nk <- SingleCellExperiment(assays = list(counts = as.matrix(do.call(cbind, nk.data))), colData = meta.nk[order(sample)])
sce.nk %<>% calculateQCMetrics()

# remove non-NK cells from data
sce.nk.filt <- subset(sce.nk,calcAverage(sce.nk) > 0)
sce.nk.filt %<>% computeSumFactors()
sce.nk.filt %<>% normalize()
assays(sce.nk.filt)[["scaled_logcounts"]] <- t(apply(assays(sce.nk.filt)[["logcounts"]],1,scale))
gnls<- c("CD3D","CD3E","CD3G","TRAC","CD8A","CD8B","GZMK","CCL5","NKG7","FCER1G","NCAM1","KLRB1","KLRD1","KLRF1","KLRC1","KLRC2","KLRC3","KLRC4","FCGR3A","FCGR3B","ITGAL","ITGAM")
set.seed(147)
sce.nk$subcluster <- factor(kmeans(t(assays(sce.nk.filt)[["scaled_logcounts"]][gnls,]),2)$cluster)

# subset and normalise data
sce.nk.filt2 <- subset(sce.nk,calcAverage(sce.nk) > 0, subcluster == 1)
sce.nk.filt2 %<>% computeSumFactors()
sce.nk.filt2 %<>% normalize()

# add scaled counts (z-score) to assays
assays(sce.nk.filt2)[["scaled_logcounts"]] <- t(apply(assays(sce.nk.filt2)[["logcounts"]],1,scale))

# identify HVGs
var.fit <- trendVar(sce.nk.filt2, use.spikes=F)
var.decomp <- decomposeVar(sce.nk.filt2, var.fit)
var.decomp$name <- rownames(var.decomp)
var.decomp <- var.decomp[with(var.decomp,order(-bio,FDR)),]
hvgs <- var.decomp[var.decomp$FDR < 0.05,'name']
## DONE