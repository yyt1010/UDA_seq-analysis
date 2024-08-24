pdf("proteimuria_Donor_enrichment.pdf",height=9,width=6)
library(ArchR)
library(Seurat)
library(ComplexHeatmap)
library(circlize)


cM=readRDS("enrich_mat_pro.rds")
cM=log(cM+0.000001,base=2)
print(max(cM))
cM1=cM
gray_scale <- colorRamp2(c(-max(cM),0,max(cM)), c("royalblue","gray","red"))

meta=read.table("meta.txt",sep="\t",header=T)
df=data.frame(donor=meta$Donor,proteinuira=meta$X24hUpro)
df$donor=paste0(df$donor,"(24hUP:",df$proteinuira,")")
df=df[df$donor %in% rownames(cM),]

df=df[order(-df$proteinuira),]

cM1=cM1[df$donor,]
print(cM1)

idx=which(substr(rownames(cM),1,2)=="WB")
cM1=rbind(cM1,cM[idx,])
print(cM1)

row_ha = rowAnnotation(proteuria = c(rep(">=3.5g",6),rep("<=3g",19),rep("Control",10))  ,col = list(proteuria = c(">=3.5g" = "brown", "<=3g" = "yellow", "Control" = "cornsilk")) )
row_ha
print(row_ha)
dim(cM1)
rownames(cM1)=substr(rownames(cM1),1,5)
p=Heatmap(cM1,heatmap_legend_param = list(
        title = "log2 \n (selected_freq/celltype_freq)"),col=gray_scale,rect_gp = gpar(col = "white", lwd = 2),cluster_rows = FALSE, cluster_columns = FALSE,right_annotation=row_ha)
print(p)

cM=readRDS("enrich_mat.rds")
print(max(cM))
p=Heatmap(cM,heatmap_legend_param = list(
        title = "enrichment_p"),col=c('red','blue'),rect_gp = gpar(col = "white", lwd = 2))
#print(p)

