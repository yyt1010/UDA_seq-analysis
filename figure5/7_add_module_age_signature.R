library(Seurat)
data=readRDS("/p300s/jiangl_group/huangzh/Post_index/PBMC/UDA1/QC_clusters/TCR/l2_new.rds")
gene=read.csv("Age_pos_Vs_None.csv")
idx=which(gene$log2FoldChange>1)
gene_table=gene[idx,]
write.table(gene_table,file="Aging_pos_gene.txt",sep="\t",row.names=F, quote=F)

gene=gene$symbol[gene$log2FoldChange>1]
features=list(gene)
data=AddModuleScore(data,features=features,name="Age_signature_pos")
saveRDS(data,file="7_add_module_pos.rds")

