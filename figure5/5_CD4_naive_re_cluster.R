library(Seurat)
data=readRDS("/p300s/jiangl_group/huangzh/Post_index/PBMC/UDA1/QC_clusters/l2_new.rds")
idx=which(data$predicted.celltype.l2 =="CD4 Naive")
data=data[,idx]
mat=data@assays$RNA@counts
saveRDS(mat,file="1_mat.rds")

mat=readRDS("1_mat.rds")
data=CreateSeuratObject(mat)

data=NormalizeData(data)
data$batch=substr(colnames(data),1,6)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(data)
data <- ScaleData(data)
print("1")
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- RunHarmony(data, group.by.vars="batch")

data <- FindNeighbors(data, reduction = "harmony", dims = 1:20)
data <- FindClusters(data, resolution = 0.5)
#data <- RunTSNE(data, dims = 1:50, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
data <- RunUMAP(data, dims = 1:20, label = T,reduction="harmony")


saveRDS(data,file="CD4_naive_harmony.rds")



