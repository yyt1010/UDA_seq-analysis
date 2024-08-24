library(Seurat)
data=readRDS("2_EC_GC_only.rds")

mat=data@assays$RNA@counts

data=CreateSeuratObject(mat)

data=NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

data <- ScaleData(data)
print("1")
data <- RunPCA(data)
data <- FindNeighbors(data,dims=1:20)
data <- FindClusters(data, resolution = 0.8)
print("2")
data <- RunUMAP(data, dims = 1:20)
saveRDS(data,file="6_after_cluster_20PC_res0.8_rep.rds")
print("3")


