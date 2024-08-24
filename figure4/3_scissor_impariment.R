t1=Sys.time()
library(Scissor)
phenotype <- readRDS("../2_impairment.rds") ##binary phenotype
phenotype  
tag <-c("non_impairment","impairment")
data <- readRDS("~/RNA/results/Scissor2/impairment/immune1/scRNA/5_after_cluster_20PC_res0.5_rep_celltype.rds")





sc_dataset=data
set.seed(123)
select=sample(1:ncol(sc_dataset),ceiling(ncol(sc_dataset)/2))
sc_dataset=sc_dataset[,select]
bulk_dataset <- readRDS("/p300s/jiangl_group/huangzh/Post_index/Kidney/Batch6/Batch6_analysis/Merge/T3/RNA/results/Scissor2/impairment/immune1/Bulk_RNA/Only_ref/1_bulk_2_count.rds")  ##bulk_profile


head(phenotype)
#sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)
DimPlot(sc_dataset, reduction = 'umap', label = T, label.size = 10)
print("1")
infos1 <- Scissor(bulk_dataset, sc_dataset,phenotype=phenotype, tag=tag, alpha = 0.05,
               family = "binomial", Save_file = 'Scissor_impairment.RData')

saveRDS(infos1$Coefs,file="3_scissor_Coefs.rds")
Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
print("2")
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
sc_dataset$Coefs <- infos1$Coefs
saveRDS(sc_dataset,file="3_scissor_Alpha0.05_scRNA_ref_count.rds")
#saveRDS(sc_dataset,file="3_scissor1.rds")
print("3")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
print("4")
t2=Sys.time()
print(t2-t1)

t2=Sys.time()
print(t2-t1)
