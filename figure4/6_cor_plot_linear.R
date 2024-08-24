library(ggplot2)
library(Seurat)
library(ggpubr)
library(ggrepel)
df=readRDS("9_df_median_immune_impairment_score.rds")
idx=which(df$New=="1")
df$donor[idx]=paste0(df$dono[idx],"_New")
#df$metabolic_impairment_score=df$immune_impairment_score
p=ggplot(df, aes(eGFR,immune_impairment_score,label=donor)) + geom_point() + geom_smooth(method="lm",se=T)+ggpubr::stat_cor(digits=3,label.x.npc = "center",
  label.y.npc = "top")+theme(panel.spacing = unit(0, "lines"))+
             theme(panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("immune_impairment")+NoLegend()+geom_text_repel(aes(label=donor),label.size=1)


#+geom_text_repel(aes(label=donor))
#+geom_label_repel(box.padding = 0.5, max.overlaps = 10)


print(p)

