library(Seurat)
atac=readRDS("/p300s/jiangl_tmp/huangzh/UDA/Species_mix_RNA_ATAC/Cellranger_human_mouse/Stat/ATAC_frag/3_proportion_atac.rds")
rna=readRDS("/p300s/jiangl_tmp/huangzh/UDA/Species_mix_RNA_ATAC/Cellranger_human_mouse/Stat/RNA/2_seurat.rds")
table(atac$barcode == colnames(rna))

idx1=which(rna$nFeature_RNA>1000)

idx_rna=idx1
idx_atac=which(atac$all>5000)

idx=unique(idx_rna,idx_atac)



atac$hg38_umi=rna$human
atac=atac[idx,]
saveRDS(atac,file="4_joint.rds")


df=atac
df$species="Multiplets"
idx1=which(df$hg38_ratio>0.85)
idx2=which(df$hg38_umi>0.85)
idx=unique(union(idx1,idx2))
df$species[idx]=paste0("human (",round(length(idx)/nrow(df),4),")")
sum=sum+length(idx)/nrow(df)



idx1=which(df$hg38_ratio < 0.15)
idx2=which(df$hg38_umi < 0.15)
idx=unique(union(idx1,idx2))
df$species[idx]=paste0("mouse (",round(length(idx)/nrow(df),4),")")
sum=sum+length(idx)/nrow(df)

df$species[df$species=="Multiplets"]=paste0("Multiplets (",round(1-sum,4),")")
p=ggplot(df,aes(x=hg38_ratio,y=hg38_umi,colour=species))+geom_point(size=0.5)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme(panel.grid=element_blank())+theme_bw()+
theme(panel.border = element_blank())+
scale_color_manual(values = c("black","red","gray"))+xlab("hg38 atac fragments ratio")+ylab("hg38 rna UMIs ratio")+
xlim(0,1)+ylim(0,1)+ggtitle("with round2 barcode")


p
p+theme(legend.position="none")

