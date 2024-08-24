library(Seurat)
library(ggplot2)
data1=readRDS("2_seurat.rds")
data=CreateSeuratObject(data1@assays$RNA@counts)
data$human=data1$human
data$mouse=data1$mouse
data$species="Multiplets"
idx_hm=which(data$human> 0.85)
idx_mm=which(data$mouse> 0.85)

data$species[idx_hm]="human"
data$species[idx_mm]="mouse"

data$Genes=data$nFeature_RNA
summary(data$Genes)
data$UMI=data$nCount_RNA
idx1=which(data$Genes>1000)
data=data[,idx1]


table(data$species)/ncol(data)
ncol(data)
saveRDS(data,file="3_filter.rds") 
