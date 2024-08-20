library(Seurat)
ref=readRDS("~/Kidney_reference/KPMP_2023_annotation/KPMP_snV3.rds")
que=readRDS("5_after_cluster_20PC_res0.5_rep.rds")

ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)


anchors <- FindTransferAnchors(reference=ref,query = que,
    dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = ref$subclass.l2
que <- AddMetaData(que, metadata = predictions)
saveRDS(que,file="4_annotated.rds")
