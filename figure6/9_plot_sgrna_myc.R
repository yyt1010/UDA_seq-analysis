library(ggplot2)
library(Matrix)
data=readRDS("1_merge_mat_mean.rds")
aa=data["MYC",]

saveRDS(aa,file="2_tmp_sgrna.rds")


df=data.frame(exp=aa)
df$label=rownames(df)
df$group="group1"
idx=which(substr(rownames(df),1,6)=="group2")
df$group[idx]="group2"

idx=which(substr(rownames(df),1,14) %in% c("group1:negCtrl","group1:Control"))
df$group[idx]="group1:Control"
idx=which(substr(rownames(df),1,14) %in% c("group2:negCtrl","group2:Control"))
df$group[idx]="group2:Control"

saveRDS(df,file="2_tmp_sgrna.rds")

ggplot(df,aes(x=group,y=exp))+geom_point()+geom_text(aes(label = label),nudge_x=0,nudge_y=0)
ggplot(df,aes(x=group,y=exp))+geom_point()
