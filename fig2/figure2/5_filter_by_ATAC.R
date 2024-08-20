##remain both high quality cell in RNA and ATAC part, then filter out donor whose cell number is less than 1000
library(Seurat)
data=readRDS("4_annotated.rds")
barcode=readRDS("~/ATAC_process/4_ATAC_QC.rds")
barcode=intersect(colnames(data),barcode)
data=data[,barcode]

aa=table(data$donor)

select=names(donor_cell[donor_cell>1000])

idx=which(data$donor %in% select)
data=data[,idx]
dim(data)
length(unique(data$donor))
saveRDS(data,file="5_add_meta_filter_1000.rds")

