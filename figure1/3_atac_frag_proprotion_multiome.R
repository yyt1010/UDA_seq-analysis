df=readRDS("atac_frag.rds")

df$hg38_reads=df$hg38_frag
df$mm10_reads=df$mm10_frag
df$all=df$hg38_reads+df$mm10_reads

df$hg38_ratio=df$hg38_frag/df$all
df$mm10_ratio=df$mm10_frag/df$all

saveRDS(df,file="3_proportion_atac.rds")
