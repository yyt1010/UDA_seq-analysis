library(Scissor)
phenotype <- readRDS("../2_proteinuria.rds")
tag <-c("non_proteinuria","proteinuria")
data <- readRDS("~/scRNA_seq/5_after_cluster_20PC_res0.5_rep_celltype.rds")

sc_dataset=data
set.seed(123)
select=sample(1:ncol(sc_dataset),ceiling(ncol(sc_dataset)/2))  #for faster calculation
sc_dataset=sc_dataset[,select]
saveRDS(sc_dataset,file="downsampling_RNA.rds")
bulk_dataset <- readRDS("~/scRNA_seq/1_bulk_mean_count.rds")


head(phenotype)
#sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)
DimPlot(sc_dataset, reduction = 'umap', label = T, label.size = 10)
print("1")
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.8,tag=tag,
                 family = "binomial", Save_file = 'Scissor.RData')

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
print("2")
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
saveRDS(sc_dataset,file="3_scissor_A_0.8_count.rds")
