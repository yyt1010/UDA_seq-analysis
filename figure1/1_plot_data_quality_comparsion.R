pdf("data_comparsion_multi_RNA.pdf")
library(ggplot2)
df=readRDS("1_df_all.rds")
idx=which(df$RNA<quantile(df$RNA,0.99)) ##0.99 quantile for a better visulization
df=df[idx,]

my_colors=c("yellow1","dodgerblue4")
p = ggplot(data=df,aes(x=sample, y=RNA)) +
  geom_violin(aes(fill=group),position="dodge",width=0.5,scale="width") + scale_fill_manual(values = my_colors)+
  geom_boxplot(outliers = F, width=0.12,outlier.colour = NA) +
  theme_bw() +
  xlab(element_blank()) +
  ylab("Gene number") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
p
ggsave("Genes.compare.png", p, width = 6, height = 8)
ggsave("PBMC_compare.pdf", p, width = 6, height = 8)


