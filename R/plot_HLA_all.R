## LOAD DATA
library(data.table)
library(Seurat)
fls <- list.files(pattern = ".h5$")
balf.data <- lapply(fls, Read10X_h5)
names(balf.data) <- sapply(strsplit(fls,"_"),"[[",2)

meta.cell <- fread("all.cell.annotation.meta.txt")[!grepl("^GSM",sample)]

# match nrows and columns
idx.all <- split(meta.cell[,gsub("_.*","",ID)],meta.cell$sample)
all.data <- mapply(balf.data[names(idx.all)], idx.all, FUN = function(x,y) x[1:min(sapply(balf.data,nrow)),gsub("-.*","",colnames(x)) %in% y, drop=F] )

library(magrittr)
library(scater)
library(scran)
sce.all <- SingleCellExperiment( assays = list(counts = as.matrix(do.call(cbind, all.data))), colData = meta.cell[order(sample)])
sce.all %<>% calculateQCMetrics()

rm(list = c("all.data","balf.data")) # cleanup

# subset and normalise data
sce.all.filt <- subset(sce.all,calcAverage(sce.all) > 0)
rm(sce.all) # cleanup

sce.all.filt %<>% computeSumFactors()
sce.all.filt %<>% normalize()

# add scaled counts (z-score) to assays
assays(sce.all.filt)[["scaled_logcounts"]] <- t(apply(assays(sce.all.filt)[["logcounts"]],1,scale))

## DONE

gnls.hla <- rownames(sce.all.filt)[grep("^HLA-",rownames(sce.all.filt))]

dat.hla <- as.data.table(melt(logcounts(sce.all.filt)[gnls.hla,]))
dat.hla[, group := rep(sce.all.filt$group,each=length(gnls.hla))]
dat.hla[, celltype := rep(sce.all.filt$celltype,each=length(gnls.hla))]

# fwrite(t(dcast(dat.hla, Var1~celltype+group+Var2, fun.aggregate=mean)), "out/hla_expression.tsv",quote = F,sep = "\t")

library(ggplot2)
library(cowplot)
p.violin.hla <- 
  ggplot(dat.hla, aes(y=value, x=group, fill=group)) +
    geom_violin(scale="width", show.legend = F) +
    labs(y="log-normalized counts") +
    facet_grid(celltype~Var1) +
    scale_fill_manual(values=c("#FFFFFF","#DE944B","#742925") ) +
    theme_cowplot(10) +
    theme(strip.background = element_blank(), axis.title.x = element_blank())

ggsave2("plots/violin_hla_all.pdf",width = 12,height = 7, p.violin.hla)

sink("out/hla_fdr.tsv")
  by(dat.hla, list(dat.hla$celltype, dat.hla$Var1), FUN = function(x) with(x, pairwise.wilcox.test(value, group, "fdr"))$p.value )
sink()
