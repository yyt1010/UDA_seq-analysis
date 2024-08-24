library(Seurat)
library(ggplot2)
library("scatterplot3d")
library(plot3D)
data=readRDS("3_filter.rds")
mat=data@assays$RNA@counts
mat=as.matrix(mat)

idx=which(substr(rownames(mat),1,6)=="GRCh38")
human_UMI=colSums(mat[idx,])

idx=which(substr(rownames(mat),1,6)=="mm10--")
mouse_UMI=colSums(mat[idx,])

idx=which(substr(rownames(mat),1,6)=="LM----")
locust_UMI=colSums(mat[idx,])


df=data.frame(human=human_UMI,mouse=mouse_UMI,locust=locust_UMI)
colors <- c("blue", "red", "yellow","grey")
names(colors)=c("human","mouse","locust","Multiplets")
data$species=factor(data$species,levels=c("human","mouse","locust","Multiplets"))
colors <- colors[data$species]
df$colors=colors
saveRDS(df,file="3_species.rds")

scatterplot3d(human_UMI,mouse_UMI,locust_UMI,angle=30,highlight.3d=TRUE)

print(as.integer(data$species))
scatter3D(human_UMI, mouse_UMI, locust_UMI, colvar = as.integer(data$species), col = c("blue", "red", "yellow","grey"), add = FALSE,colkey = list(at = c(1,2, 3, 4), side = 1, 
          addlines = TRUE, length = 0.5, width = 0.5,labels = c("human (0.3752)","mouse(0.2805)","locust(0.3232)","Multiplets(0.0211)")),pch = 16,theta = 120, main = "Species mixture", xlab = "Human UMIs",
          ylab ="Mouse UMIs", zlab = "Locust UMIs")
