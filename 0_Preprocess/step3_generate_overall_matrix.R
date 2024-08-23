library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
length=96  ##index number, typically 96. 
datalist=list()
for (i in 1:length)
{  
    path=paste0("~/cellranger_results_path/Index_",i,"/Index_",i,"/outs/filtered_feature_bc_matrix")  ##read matrix for each index
    print(path)
    datalist[[i]]= Read10X(data.dir = path)$Gene ##if not multi-ome, datalist[[i]]= Read10X(data.dir = path)
    datalist[[i]]=CreateSeuratObject(counts = datalist[[i]], min.cells = 0, min.features = 0) ##no filtering here
    
   
}
data=c(datalist[[1]])
for (i in 2:length)
{
    data[i]=datalist[[i]]
}
merged_data=merge(x=datalist[[1]],y=data[2:length],add.cell.ids=1:96)


##add demuxlet results

barcode=c()
donor=c()
SNG=c()
SNP_num=c()
SNG_POSTERIOR=c()

for (i in 1:length){
    path=paste0("~/path_demuxlet_result/Index_",i,"/",i,"_result_GP.best")}
    deplex=read.table(path,header=T)
    barcode=c(barcode,paste0(i,"_",deplex$BARCODE))
    donor=c(donor,deplex$SNG.BEST.GUESS)
    SNG=c(SNG,deplex$DROPLET.TYPE)
    SNP_num=c(SNP_num,deplex$NUM.SNPS)
    SNG_POSTERIOR=c(SNG_POSTERIOR,deplex$SNG.POSTERIOR)
}
merged_data$donor=donor
merged_data$SNG=SNG
merged_data$SNP_num=SNP_num
merged_data$POSTERIOR=merged_data$POSTERIOR



##filter low quality cells

merged_data[["percent.mt"]] <- PercentageFeatureSet(merged_data, pattern = "^MT-")
merged_data@meta.data$UMI=merged_data@meta.data$nCount_RNA
merged_data@meta.data$Genes=merged_data@meta.data$nFeature_RNA


idx=which(merged_data$Genes > 200)   ##usr can set suitable  cutoff 
merged_data=merged_data[,idx]
merged_data <- subset(merged_data, percent.mt < 5)


saveRDS(merged_data,file="merged_gene_mat.rds")

