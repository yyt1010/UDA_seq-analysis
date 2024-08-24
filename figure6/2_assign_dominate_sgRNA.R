library(Matrix)
data=readRDS("domi_mat.rds")
Max=apply(data,2,max)
idx=which(Max>0.51)
print(length(idx))
data=data[,idx]
sgRNA=vector(mode="numeric",length=ncol(data))
for (i in 1:length(sgRNA)){
    sgRNA[i]=rownames(data)[which(data[,i]>0.51)]
}
names(sgRNA)=colnames(data)
saveRDS(sgRNA,file="2_assign.rds")
