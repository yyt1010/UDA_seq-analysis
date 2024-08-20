library(Seurat)

##QC by gene number, MT-ratio, SNP number and singlets. 
data=readRDS("RNA_seurat.rds")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data@meta.data$UMI=data@meta.data$nCount_RNA
data@meta.data$Genes=data@meta.data$nFeature_RNA
idx=which(data$Genes > 200)
data=data[,idx]
data <- subset(data, percent.mt < 5)
idx=which(data$Drop=="SNG")
data=data[,idx]
idx=which(data$SNP_num>50)
data=data[,idx]



##clustering and umap embedding

data=CreateSeuratObject(data)
data=NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

data <- ScaleData(data)

data <- RunPCA(data)
data <- FindNeighbors(data,dims=1:20)
data <- FindClusters(data, resolution = 0.5)

data <- RunUMAP(data, dims = 1:20)
saveRDS(data,file="5_after_cluster_20PC_res0.5.rds")


