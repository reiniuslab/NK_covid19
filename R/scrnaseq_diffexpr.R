source("R/load_data_scran_NK.R")

##

pairwise.zlm <- function(sce, group, threshold = 0.2){
  require(MAST)
  require(data.table)
  sca <- SceToSingleCellAssay(sce)
  v <- combn(unique(group),2,simplify = F)
  
  # performs pairwise differential expression on contrasts based on {group}
  ls <- lapply(v,function(x){
    # subset data
    idx <- which(group %in% x)
    expressed <- freq(sca[,idx]) > threshold
    
    sca.red <- sca[expressed,idx]
    colData(sca.red)$cngeneson <- scale(colSums(assay(sca.red)>0))
    colData(sca.red)$group <- factor(group[idx])
    
    # perform differential expression
    fit <- zlm(~group + cngeneson, sca.red)
    hypothesis <- paste0("group",x[2])
    res <- summary(fit,doLRT=hypothesis )$datatable
    
    # format output
    out <- merge(res[contrast==hypothesis & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 res[contrast==hypothesis & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    out[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    setnames(out, c("primerid","coef","Pr(>Chisq)"), c("gene","logFC","p-value"))
    setorder(out,fdr)
    return(out)
  })
  names(ls) <- lapply(v,function(x) paste0(rev(x),collapse="_v_") )
  return(ls)
}

##

res.ls <- pairwise.zlm(sce.nk.filt2,sce.nk.filt2$group)
res.dt <- rbindlist(res.ls, idcol=T)

# write output
lapply(names(res.ls), function(x) fwrite(res.ls[[x]],paste0("out/",x,".tsv"), quote = F, sep = "\t") )

gn.diff <- res.dt[fdr < 1e-3, unique(gene[order(fdr)])]

# factoextra::fviz_nbclust(assay(sce.nk.filt2,"scaled_logcounts")[res.dt[fdr < 1e-3, unique(gene)],], kmeans, "gap_stat") # 6 clusters
set.seed(42) # reproducible cluster order
gn.cl <- kmeans(assay(sce.nk.filt2,"scaled_logcounts")[gn.diff,], 6, nstart=25, iter.max = 100)$cluster
gn.topcl <- unlist(tapply(gn.diff,gn.cl,function(x) head(x,length(x)*0.1) ))

write.table(sort(gn.cl), file = "out/gene_clusters.tsv", quote = F,col.names = F,sep = "\t")

dat.melt <- as.data.table(melt(assay(sce.nk.filt2,"scaled_logcounts")[gn.diff,]))
dat.melt[,cl := rep(gn.cl,ncol(sce.nk.filt2))]
dat.melt[,group := rep(sce.nk.filt2$group, each = length(gn.diff))]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
p.heat.diffexprs <-
  Heatmap(
    matrix = assay(sce.nk.filt2,"scaled_logcounts")[gn.diff,],
    name="Z-score",
    show_column_names = F,
    show_row_names = F,
    show_row_dend = F,
    show_column_dend = F,
    column_split = sce.nk.filt2$group,
    row_split = gn.cl,
    col=colorRamp2(seq(-4,4,length.out = 11), rev(brewer.pal(11,"PiYG"))),
    right_annotation = rowAnnotation(foo = anno_mark(at = match(gn.topcl,row.names(assay(sce.nk.filt2,"scaled_logcounts")[gn.diff,])), labels = gn.topcl))
  )

pdf("plots/heatmap_differential_expression.pdf")
  draw(p.heat.diffexprs)
dev.off()

library(ggplot2)
library(cowplot)
p.linepoint.clust <-
  ggplot(dat.melt[,mean(value,na.rm=T),by=c("Var1","cl","group")], aes(y=V1, x=group)) +
  stat_summary(fun.y="mean", geom="line", aes(group=cl)) +
  stat_summary(fun.data="mean_cl_normal") +
  facet_wrap(~cl, nrow=1) +
  labs(y="Average Z-score", x=NULL) +
  theme_cowplot(10) +
  theme(strip.background = element_blank())

ggsave2("plots/linepoint_cluster.pdf",height = 2, width = 7, p.linepoint.clust)

## 
library(ggplot2)
library(cowplot)
dat.melt$logcounts <- as.data.table(melt(assay(sce.nk.filt2,"logcounts")[gn.diff,]))$value

dat.melt.agg <- dat.melt[,list(z = mean(value), det = mean(logcounts > 0)), by=c("cl", "group")]
dat.melt.agg[,z_fix := z]
dat.melt.agg[z_fix > 0.5,z_fix := 0.5]
dat.melt.agg[z_fix < -0.5,z_fix := -0.5]

p.point.clust <- 
  ggplot(dat.melt.agg, aes(x=factor(cl), y=group, col=z_fix, size=det)) +
  geom_point() +
  coord_flip() +
  labs(col="Z-score",size="Detection") +
  scale_x_discrete(labels=c("(1) Effector function/\nchemotaxis", "(2) Activation/\nproliferation", "(3) Metabolism/\ndegranulation", "(4) Interferon\nresponse","(5) Mitochondrial\nactivity", "(6) Translation/\nRNA metabolism")) +
  scale_colour_gradient2(low = "#A5A4A4", mid = "#FFFFFF", high = "#683B97") +
  #scale_color_distiller(palette="RdBu", limits=c(-0.5,0.5) ) +
  theme_cowplot(10) +
  theme(axis.title = element_blank())

ggsave2("plots/pointheat_cluster.pdf",height = 3, width = 3, p.point.clust)

## Enrichment for signature genes
library(MAST)
sca.nk <- SceToSingleCellAssay(sce.nk.filt2)
sca.nk.filt <- sca.nk[freq(sca.nk) > 0.2]

colData(sca.nk.filt)$cngeneson <- scale(colSums(assay(sca.nk.filt)>0))

# perform differential expression
fit.nk <- zlm(~group + cngeneson, sca.nk.filt)

# boots.nk <- bootVcov1(fit.nk, 100)
# saveRDS(boots.nk, "out/boots.nk.rds")
boots.nk <- readRDS("out/boots.nk.rds")

gnls.signature <- fread("41467_2019_11947_MOESM6_ESM.tsv")
set <- with(gnls.signature, split(gene,paste0("cl",cluster) ))
set.idx <- limma::ids2indices(set, rownames(sca.nk.filt) )

gseaM <- gseaAfterBoot(fit.nk, boots.nk, set.idx, CoefficientHypothesis("groupM"))
gseaS <- gseaAfterBoot(fit.nk, boots.nk, set.idx, CoefficientHypothesis("groupS"))

z.dt <- rbindlist(idcol = T, list(
  "M" = as.data.table(calcZ(gseaM,combined = "stouffer"),keep.rownames = T),
  "S" = as.data.table(calcZ(gseaS,combined = "stouffer"),keep.rownames = T)
))

p.gsea.nk <- 
  ggplot(z.dt, aes(x=rn, y=Z, size=-log10(P), col = .id )) +
    geom_point() +
    labs(y="Z-score", x=NULL, col=NULL) +
    coord_cartesian(ylim=c(-8,8)) +
    theme_cowplot(10) +
    scale_color_manual(values=c("#DE944B","#742925"))

ggsave2("plots/gsea_nk_signatures.pdf", width = 3,height = 3, p.gsea.nk)
