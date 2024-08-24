library(Seurat)
data=readRDS("after_cluster_20PC_res0.5_harmony.rds")
DimPlot(data,raster=F,label=T,repel=T)+NoLegend()
idx=which(!(data$seurat_clusters %in% c("4","9","12"))) ##low quality cell clusters
data=data[,idx]
DimPlot(data,raster=F,label=T,repel=T)+NoLegend()

dim(data)
mat=data@assays$RNA@counts
saveRDS(mat,file="2_mat.rds")

