library(Seurat)

data=readRDS("merged_gene_mat.rds")
doublet=read.table("scru_doublet.txt")
idx=which(doublet=="False")
data$DROP_type_scr="DBL"
data$DROP_type_scr[idx]="SNG"


idx1=which(data$SNP_num >50)
idx2=which(data$SNG_POSTERIOR>0.99)
idx=intersect(idx1,idx2)
data=data[,idx]
data$final_drop="singlet"

idx1=which(data$type=="DBL")
idx2=which(data$DROP_type_scr=="DBL")
idx_DBL=intersect(idx1,idx2)

data$final_drop[idx_DBL]="doublet"

idx=which(data$final_drop=="singlet")
data=data[,idx]
saveRDS(data,file="5_seurat_final.rds")

