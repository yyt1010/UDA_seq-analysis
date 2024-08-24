library(Seurat)
library(ggplot2)
data=readRDS("6_add_module_impairment1_Alpha0.05_scRNA.rds")
data$impairment_signature=data$proteinuria1
meta=read.table("meta.txt",sep="\t",header=T)




FeaturePlot(data,features="impariment_signature",cols=c("blue","red"),raster=T)
bb=table(data$predicted.id)
select=names(bb[bb>200])
idx=which(data$predicted.id %in% select)
data=data[,idx]
impairment_signature=data$impariment_signature
donor=data$donor
celltype=data$predicted.id



df=data.frame(impairment=impariment_signature,donor=donor,celltype=celltype)


VlnPlot(data,features="impairment_signature",group.by="predicted.id",pt.size=0,sort=T)+NoLegend()+geom_boxplot(width=0.3, fill="white")
VlnPlot(data,features="impairment_signature",group.by="donor",pt.size=0,sort=T)+NoLegend()+geom_boxplot(width=0.3, fill="white")


data$donor=factor(data$donor,levels=donor_order)
data$celltype=factor(data$predicted.id,levels=celltype_order)
p1=VlnPlot(data,features="impairment_signature",group.by="celltype",pt.size=0,sort=F)+NoLegend()+geom_boxplot(width=0.3, fill="white",outlier.shape=NA)
p2=VlnPlot(data,features="impairment_signature",group.by="donor",pt.size=0,sort=F)+NoLegend()+geom_boxplot(width=0.3, fill="white",outlier.shape=NA)
p1+ggtitle("impariment_signature")
p2+ggtitle("impariment_signature")

