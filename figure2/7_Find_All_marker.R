pdf("markers_Kidney.pdf",height=20,width=40)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
data=readRDS("5_add_meta_filter_1000.rds")
idx=which(!is.na(data$celltype_l1))
data=data[,idx]
Idents(object=data)=data$celltype_l1
if (TRUE){
data.markers <- FindAllMarkers(data,only.pos = TRUE,logfc.threshold = 0.25,max.cells.per.ident=1000)
saveRDS(data.markers,file="7_all_markers.rds")

gene=readRDS("7_all_markers.rds")
#idx=which(gene$cluster %in% celltype)
#gene=gene[idx,]
#gene$cluster=factor(gene$cluster,levels=celltype)

gene %>%
    group_by(cluster) %>%

    top_n(n = 2, wt = avg_log2FC) -> top10

saveRDS(top10,file="top_gene.rds")
}
top10=readRDS("top_gene.rds")
data=ScaleData(data,features=top10$gene)
#aa=table(data$predicted.id)
#select=names(aa[aa>1000])
#data=data[,data$predicted.id %in% celltype]

DoHeatmap(subset(data,downsample=100),features = top10$gene) + NoLegend()
DoHeatmap(subset(data,downsample=100),features = top10$gene)
aa=make.unique(top10$gene)

levels(data)=celltype
DotPlot(data,features=aa)+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

