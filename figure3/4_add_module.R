library(Seurat)
data=readRDS("../../../10_add_meta_flter_1000.rds")
gene=readRDS("5_gene_no_gender.rds")
features=list(gene)
data=AddModuleScore(data,features=features,name="proteinuria")




saveRDS(data,file="4_add_module_proteinuria1_A_0.8_count.rds")
