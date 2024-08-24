library(rtracklayer)
data=readRDS("2_markers.rds")
idx1=which(data$p_val_adj<0.01)
idx2=which(data$avg_log2FC >1)
idx=intersect(idx1,idx2)

gene=rownames(data)[idx]


gtf=import("/xtdisk/jiangl_group/huangzh/software/software/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes1.gtf")
idx=which(gtf@seqnames =="chrY")
gene_chrY=unique(gtf$gene_name[idx])

gtf=import("/xtdisk/jiangl_group/huangzh/software/software/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes1.gtf")
idx=which(gtf@seqnames =="chrX")
gene_chrX=unique(gtf$gene_name[idx])

gender=c(gene_chrY,gene_chrY)

gene=gene[!(gene %in% gender)]


idx=which(rownames(data) %in% gene)
data=data[idx,]
write.table(data,file="5_proteinuria_pos.txt",sep="\t",quote=F)
saveRDS(gene,file="5_gene_no_gender.rds")
