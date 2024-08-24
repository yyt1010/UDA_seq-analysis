library(ArchR)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

data=readRDS("/p300s/jiangl_group/huangzh/Post_index/Kidney/Batch6/Batch6_analysis/Merge/T3/RNA/results/Scissor1/2_proteinuria_3g_rep/All_cell_0.8/3_scissor_A_0.8_count.rds")


cM <- confusionMatrix(paste0(data$donor), paste0(data$scissor))
cM=as.matrix(cM)
colnames(cM)=c("Backgroud","Proteinuria_pos","Proteinuria_neg")
saveRDS(cM,file="1_donor_mat.rds")

bb=colSums(cM)
for (i in 1:ncol(cM)){
    cM[,i]=cM[,i]/bb[i]
}
saveRDS(cM,file="1_donor_mat_pro.rds")
