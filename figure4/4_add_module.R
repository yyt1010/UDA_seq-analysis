library(Seurat)
data=readRDS("..//scRNA/1_I_Control.rds")

gene=readRDS("5_gene_no_gender.rds")




features=list(gene)
data=AddModuleScore(data,features=features,name=impairment")




saveRDS(data,file="6_add_module_impairment1_Alpha0.05_scRNA.rds")


