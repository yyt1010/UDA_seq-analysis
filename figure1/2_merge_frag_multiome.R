df=readRDS("/p300s/jiangl_tmp/huangzh/UDA/Species_mix_RNA_ATAC/Cellranger_human_mouse/Index_1/index_1/outs/1_frag_counts.rds")
df$barcode=paste0("post_index",1,"_",df$barcode)
for (i in 2:16){
print(i)
path=paste0("/p300s/jiangl_tmp/huangzh/UDA/Species_mix_RNA_ATAC/Cellranger_human_mouse/Index_",i,"/index_",i,"/outs/1_frag_counts.rds")
data=readRDS(path)
data$barcode=paste0("post_index",i,"_",data$barcode)
df=rbind(df,data)
}
saveRDS(df,file="atac_frag.rds")
