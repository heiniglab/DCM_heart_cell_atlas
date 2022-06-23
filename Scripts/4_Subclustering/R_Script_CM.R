library(Seurat)
library(Matrix)
library(SeuratDisk)
library(harmony)
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)


#load overall object#
#subset cardiomyocytes#

mito.features <- grep(pattern = "^MT-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = pbmc, slot ='counts'))
pbmc[['percent_mito']] <- percent.mito

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress=c("percent_mito", "n_counts"))
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE, npcs=50)
ElbowPlot(object = pbmc, ndims = 50)

options(repr.plot.height = 2.5, repr.plot.width = 6)
pbmc <- pbmc %>% RunHarmony(c("orig.ident"), plot_convergence = TRUE, epsilon.harmony=-Inf,  epsilon.cluster=-Inf)

pbmc <- pbmc %>%
  RunUMAP( reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50,  n.trees = 50, k.param = 15) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
  
DimPlot(pbmc, reduction = "umap", pt.size = .1)
ggsave("UMAP_harmony_hvg3000_pc50_ntree50_kparam15_res0.8.png")

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("Heatmap_harmony_hvg3000_pc50_ntree50_kparam15_res0.8.pdf",  width = 7, height = 15)

VlnPlot(pbmc, features = c("n_genes", "n_counts", "percent_mito", "percent_ribo", "solo_score"), ncol = 5)
ggsave("QC_harmony_hvg3000_pc50_ntree50_kparam15_res0.8.png")

saveRDS(pbmc,file="Heatmap_harmony_hvg3000_pc50_ntree50_kparam15_res0.8.rds")

#annotation#
Cluster0 = vCM2
Cluster1 = vCM1.0
Cluster2 = vCM3.0
Cluster3 = vCM1.0
Cluster4 = vCM1.0
Cluster5 = vCM1.1
Cluster6 = vCM1.2
Cluster7 = vCM1.0
Cluster8 = vCM1.3
Cluster9 = vCM4
Cluster10 = vCM3.1
Cluster11 = vCM1.3
Cluster12 = vCM5
Cluster13 = vCM1.3
Cluster14 = vCM1.0
Cluster15 = vCM3.0
Cluster16 = vCM1.0

saveRDS(pbmc,file="CM_hvg3000_pc50_ntree50_kparam15_annotated.rds")

SaveH5Seurat(pbmc3k.final, filename = "CM_hvg3000_pc50_ntree50_kparam15_annotated.h5Seurat")
Convert("CM_hvg3000_pc50_ntree50_kparam15_annotated.h5Seurat", dest = "h5ad")




