source("R/load_data_scran_NK.R")

## DIMRED
set.seed(192) # reproducibility
sce.nk.filt2 %<>% runUMAP(feature_set = head(hvgs,1000))

library(cowplot)
ggsave2("plots/dim_reduction.pdf", width = 4.8,height = 4,
  ggplot(as.data.frame(reducedDim(sce.nk.filt2, "UMAP")), aes(x=V1, y=V2, fill=sce.nk.filt2$group)) +
    geom_point(shape = 21, size=3)  +
    labs(x="UMAP 1",y="UMAP 2")  +
    scale_fill_manual(name = "Group",values=c("#FFFFFF","#DE944B","#742925")) +
    theme_cowplot(10)
)

ggsave2("plots/dim_reduction_ind.pdf", width = 4.8,height = 4,
  ggplot(as.data.frame(reducedDim(sce.nk.filt2, "UMAP")), aes(x=V1, y=V2, fill=sce.nk.filt2$sample_new)) +
    geom_point(shape = 21, size=3)  +
    labs(x="UMAP 1",y="UMAP 2")  +
    scale_fill_manual(name="Sample", values=c("#67bd63","#5b3788","#b4ad3d","#6971d7","#71883b","#c972c4","#45c097","#b1457b","#c17d38","#6d8dd7","#b74e37","#ba4758")) +
    theme_cowplot(10)
)