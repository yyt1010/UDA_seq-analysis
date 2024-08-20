library(Seurat)
pdf("Kidney_batch6_basic.pdf")
library(ggplot2)
library(Seurat)
data=readRDS("5_add_meta_filter_1000.rds")
DimPlot(data,label=T,repel=T,raster=T,shuffle=T,group.by="predicted.id")+NoLegend()+ggtitle(paste0("Celltype: ",ncol(data)," cells"))

