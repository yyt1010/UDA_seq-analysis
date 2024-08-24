pdf("4_advanced_clonal.pdf")
###colors
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


library(scRepertoire)
library(Seurat)
library(ggplot2)
if (TRUE){
combined=readRDS("1_combined_rep.rds") 

for (i in 1:nrow(combined[[1]])){
    print(i)
    combined[[1]]$barcode[i]=substr(combined[[1]]$barcode[i],15,nchar(combined[[1]]$barcode[i]))
}

print("1")
seurat=readRDS("integration4_1_50_add_donor.rds")
idx=which(seurat$predicted.celltype.l2.score>0.4)
seurat=seurat[,idx]


seurat <- combineExpression(combined, seurat, 
            cloneCall="gene+nt", group.by = "sample")
print("2")
saveRDS(seurat,file="4_seurat.rds")
}

seurat=readRDS("4_seurat.rds")
seurat@meta.data$CellType=seurat@meta.data$predicted.id
DimPlot(seurat, group.by = "cloneType",raster=F) +
    scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())
clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "donor") + 
                guides(color = "none")
print("HZ1")
library(ggraph)
#No Identity filter
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "predicted.id",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "gene+nt")
occupiedscRepertoire(seurat, x.axis = "CellType")+theme(axis.text.x=element_text(angle=45,hjust=1))
occupiedscRepertoire(seurat, x.axis = "donor")+theme(axis.text.x=element_text(angle=45,hjust=1))

alluvialClonotypes(seurat, cloneCall = "gene+nt", 
                   y.axes = c("donor","CellType"),color="CellType")

