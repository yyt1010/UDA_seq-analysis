library(Seurat)
library(Matrix)

data=readRDS("1_add_group.rds")
str(data@meta.data)
#data$gene=data1$gene
#data$sgRNA=data1$sgRNA

idx=which(data$gene %in% c("negCtrl","Control"))
data$gene[idx]="Control"




#idx=which(data$gene %in% rownames(data))
#data=data[,idx]
#data$gene=paste0(data$group,":",data$gene)
sgRNA=unique(data$sgRNA)

print(length(sgRNA))
merge_mat=matrix(data = NA, nrow = nrow(data), ncol = length(sgRNA), byrow = FALSE,dimnames = NULL)
for (i in 1:length(sgRNA)){
    idx=which(data$sgRNA==sgRNA[i])
    mat=data@assays$RNA@counts[,idx]
    #merge_mat[,i]=apply(mat,1,median)
    merge_mat[,i]=apply(mat,1,mean)
}
rownames(merge_mat)=rownames(data)
colnames(merge_mat)=sgRNA
saveRDS(merge_mat,file="1_merge_mat_mean_sgRNA.rds")
