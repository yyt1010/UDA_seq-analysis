##groups were based on pca space

library(Seurat)
data=readRDS("../../4_add_target_gene.rds")
group1=c("8","0","7")
group2=c("1","2","3","4","5","6")
idx=which(data$seurat_clusters %in% c(group1,group2))
data=data[,idx]

data$group="group1"
idx=which(data$seurat_clusters %in% group2)
data$group[idx]="group2"
saveRDS(data,file="1_add_group.rds")
