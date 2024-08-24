pdf("MYC_heatmap_down1.pdf",width=15,height=10)

library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(Matrix)
data=readRDS("6_mat_down.rds")
print(colnames(data))
data1=data[1:48,setdiff(1:23,19)]
rownames(data1)=gsub("group1:","",rownames(data1))
#colnames(data1)=paste0("1_",colnames(data1))

data2=data[49:96,setdiff(1:23,19)]
rownames(data2)=gsub("group1:","",rownames(data2))
#colnames(data2)=paste0("2_",colnames(data2))


mat=cbind(data1,data2)

p=Heatmap(mat,na_col = "black",cluster_rows = T,cluster_columns=F,heatmap_legend_param = list(
        title = "Enrichment_score"),row_title="Target gene",rect_gp = gpar(col = "white", lwd = 1),column_split = rep(c("group1", "group2"), each=22),name="mat",col=c("lightgray","red","red4"))

p

