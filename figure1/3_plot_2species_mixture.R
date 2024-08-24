library(Seurat)
library(ggplot2)
library(Matrix)
data=readRDS("3_filter.rds")
count_mat=data@assays$RNA@counts
idx_hm_gene=which(substr(rownames(data),1,2)=="GR")
idx_mm_gene=which(substr(rownames(data),1,2)=="mm")
hm_counts=colSums(count_mat[idx_hm_gene,])
mm_counts=colSums(count_mat[idx_mm_gene,])
species=vector(mode="character",length=ncol(data))

idx_hm=which(data$species=="human")
idx_mm=which(data$species=="mouse")
species[1:length(species)]=paste0("Multiplets (",round(1-(length(idx_hm)+length(idx_mm))/length(species),4),")")
species[idx_hm]=paste0("human (",round(length(idx_hm)/length(species),4),")")
species[idx_mm]=paste0("mouse (",round(length(idx_mm)/length(species),4),")")



df=data.frame(human=hm_counts,mouse=mm_counts,type=species)
p=ggplot(df,aes(x=human,y=mouse,colour=type))+geom_point(size=0.5)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme(panel.grid=element_blank())+theme_bw()+theme(panel.border = element_blank())+
scale_color_manual(values = c("black","red","gray"))+xlab("hg38 UMIs per cell")+ylab("mm10 UMIs per cell")+
xlim(0,max(c(df$human,df$mouse)))+ylim(0,max(c(df$human,df$mouse)))+ggtitle("with round2 barcode")

p
p+theme(legend.position="none")
