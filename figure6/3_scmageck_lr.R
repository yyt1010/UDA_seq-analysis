library(scMAGeCK)
library(Seurat)
# set the BARCODE and RDS file path 
BARCODE = "./3_df.txt"  ##sgRNA and cell barcode dataframe
## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = "../after_cluster_20PC_res0.5.rds"  ##scRNA_seq

# Run scmageck_lr function
# By default, the result will be saved to the current working directory. 
lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr_new1', 
                         NEGCTRL = c('Control'), PERMUTATION = 100)
lr_score <- lr_result[[1]]
lr_score_pval <- lr_result[[2]]

saveRDS(lr_result,file="4_lr_results.rds")
