library(Seurat)
data=readRDS("2_wnn_pre.rds") ##add ArchR reduction in seurat_obj


data <- FindMultiModalNeighbors(
  data, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight"
)

data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


saveRDS(data,file="8_wnn.rds")
p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
p2 <- DimPlot(data, reduction = 'wnn.umap', group.by = 'predicted.id', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
p1 
p2
