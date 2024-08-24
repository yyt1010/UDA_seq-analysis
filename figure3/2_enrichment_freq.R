
library(Matrix)
mat=readRDS("1_donor_mat.rds")
all=sum(mat)
celltype=rowSums(mat)
cell_pro=celltype/all
table(cell_pro)
data=readRDS("1_donor_mat_pro.rds")
print(data)
mat_enrich=data
for (i in 1:nrow(data)){
    for (j in 1:ncol(data)){

        
        mat_enrich[i,j]=data[i,j]/cell_pro[rownames(data)[i]]
    }



}
print(mat_enrich)
#mat_enrich=-log(mat_enrich,10)
print(mat_enrich)
saveRDS(mat_enrich,file="enrich_mat_pro.rds")
