library(clusterProfiler)
library(org.Hs.eg.db)

data=readRDS("../1_add_group.rds")
idx=which(data$gene %in% c("negCtrl","Control"))
data$gene[idx]="Control"

data$gene=paste0(data$group,":",data$gene)
if (FALSE){
idx=which(data$group=="group1")
data1=data[,idx]
tgene=setdiff(unique(data1$gene),"group1:Control")
for (i in 1:length(tgene)){
    print(i)
    path1=paste0("./Genes_down/",tgene[i],".rds")
    gene=readRDS(path1)
    tmp=select(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL')$ENTREZID
    ck=enrichKEGG(tmp)
    path2=paste0("./Kegg_down/",tgene[i],".rds")
    saveRDS(ck,file=path2)

}
}

idx=which(data$group=="group2")
data1=data[,idx]
tgene=setdiff(unique(data1$gene),"group2:Control")
for (i in 1:length(tgene)){
    print(i)
    path1=paste0("./Genes_down/",tgene[i],".rds")
    gene=readRDS(path1)
    tmp=select(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL')$ENTREZID
    ck=enrichKEGG(tmp)
    path2=paste0("./Kegg_down/",tgene[i],".rds")
    saveRDS(ck,file=path2)
}


