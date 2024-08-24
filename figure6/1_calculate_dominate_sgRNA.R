library(Seurat)
library(Matrix)
gRNA=readRDS("../../4_gRNA.rds")
data=readRDS("../T2/after_cluster_20PC_res0.5.rds")

gRNA=gRNA[,colnames(data)]

mat=gRNA@assays$RNA@counts
mat=mat+0.001
str(mat)
mat=as.matrix(mat)
for (i in 1:ncol(mat)){
    mat[,i]=(mat[,i]/sum(mat[,i]))
}
saveRDS(mat,file="domi_mat.rds")
