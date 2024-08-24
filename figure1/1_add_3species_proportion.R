library(Seurat)
library(Matrix)
data=readRDS("merged_gene_mat.rds")
idx_hm_gene=which(substr(rownames(data),1,2)=="GR")
idx_mm_gene=which(substr(rownames(data),1,2)=="mm")
idx_lm_gene=which(substr(rownames(data),1,2)=="LM")

count_mat=data@assays$RNA@counts
All_counts=colSums(count_mat)

hm_counts=colSums(count_mat[idx_hm_gene,])
mm_counts=colSums(count_mat[idx_mm_gene,])
lm_counts=colSums(count_mat[idx_lm_gene,])

data$human=hm_counts/All_counts
data$mouse=mm_counts/All_counts
data$locust=lm_counts/All_counts
saveRDS(data,file="2_seurat.rds")
