##/software/biosoft/software/python/anaconda3-python3-2018/bin/R
library(Seurat)
library(patchwork)
library(harmony)
data=readRDS("2_mat.rds")
data=CreateSeuratObject(data)
#idx=which(data$nFeature_RNA>400)
#data=data[,idx]
#dim(data)
data$batch=substr(colnames(data),1,6)
data=NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(data)
data <- ScaleData(data)
print("1")
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- RunHarmony(data, group.by.vars="batch")
saveRDS(data,file="after_cluster_20PC_res0.5_harmony.rds")
data <- FindNeighbors(data, reduction = "harmony", dims = 1:20)
data <- FindClusters(data, resolution = 0.5)
#data <- RunTSNE(data, dims = 1:50, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
data <- RunUMAP(data, dims = 1:50, label = T,reduction="harmony")
saveRDS(data,file="after_cluster_20PC_res0.5_harmony1.rds")

##annotation
reference=LoadH5Seurat("/xtdisk/jiangl_group/huangzh/SuperLoading/T20/Cellranger_intron/Stat/Seurat/Intergation/pbmc_multimodal.h5seurat")
query=readRDS("after_cluster_20PC_res0.5_harmony1.rds")
#query <- SCTransform(query, verbose = FALSE)
if (TRUE){
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reduction="pcaproject",
  reference.reduction="spca",
  dims = 1:50
)
saveRDS(anchors,file="anchors.rds")
}
anchors=readRDS("anchors.rds")
query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"),
  reference.reduction="spca",
  reduction.model = "wnn.umap"
)
saveRDS(query,file="integration4_1_50.rds")

