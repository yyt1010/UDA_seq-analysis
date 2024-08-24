for i in {3..16} #1..16
do
echo ${i}
cd /p300s/jiangl_tmp/huangzh/UDA/Species_mix_RNA_ATAC/Cellranger_human_mouse/Index_${i}/index_${i}/outs
gunzip atac_fragments.tsv.gz
cat > frag_count.R << EOF
barcode=read.table("./filtered_feature_bc_matrix/barcodes.tsv")\$V1
data=read.table("atac_fragments.tsv")
idx=which(data\$V4 %in% barcode)
data=data[idx,]

hg38=vector(mode="numeric",length=length(barcode))
mm10=vector(mode="numeric",length=length(barcode))

for (i in 1:length(barcode)){
    print(i)
    idx=which(data\$V4==barcode[i])
    tmp=data[idx,]
    idx_hm=which(substr(tmp\$V1,1,4)=="hg38")
    hg38[i]=sum(tmp\$V5[idx_hm])
    idx_mm=which(substr(tmp\$V1,1,4)=="mm10")
    mm10[i]=sum(tmp\$V5[idx_mm])
}

df=data.frame(barcode=barcode,hg38_frag=hg38,mm10_frag=mm10)
saveRDS(df,file="1_frag_counts.rds")


EOF
Rscript frag_count.R

done
