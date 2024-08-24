pdf("5_advanced_clonal.pdf")
###colors
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF",
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


library(scRepertoire)
library(Seurat)
library(ggplot2)
seurat=readRDS("l2_new.rds")



str(seurat@meta.data)
library(circlize)
library(scales)

seurat$donor=paste0(seurat$age,"_",seurat$donor)


table(seurat$donor)
p=clonalDiversity(seurat, cloneCall = "gene+nt", n.boots=100,group.by="donor")
saveRDS(p,file="clonalDiversity.rds")
print("OK")
clonalHomeostasis(combined2, cloneCall = "gene+nt")+theme(axis.text.x=element_text(angle=45,hjust=1))
clonalProportion(combined2, cloneCall = "gene+nt")
clonalOverlap(combined2, cloneCall="gene+nt", method="overlap")+theme(axis.text.x=element_text(angle=45,hjust=1))
occupiedscRepertoire(seurat, x.axis = "range",proportion=T)

vizGenes(combined2, gene = "V",
         chain = "TRB",
         plot = "heatmap",
         order = "variance",
         scale = TRUE)
vizGenes(combined2, gene = "V",
         chain = "TRA",
         plot = "heatmap",
         order = "variance",
         scale = TRUE)

clonesizeDistribution(combined2,
                      cloneCall = "gene+nt",
                      method="ward.D2")
occupiedscRepertoire(seurat, x.axis = "range",proportion=T)


