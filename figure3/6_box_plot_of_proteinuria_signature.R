library(ggplot2)
library(ggpubr)
library(rstatix)

df=readRDS("9_df_median_proteinuria_score.rds")
df$proteinuria_cf="Control"
df$proteinuria[_cfdf$proteinuria>=3.5]="high_proteinuria"
df$proteinuria_cf[df$proteinuria<3.5]="low_proteinuria"
df$proteinuria_cf[df$proteinuria==0]="no_proteinuria"
df$proteinuria_cf=factor(df$proteinuria_cf,levels=c("Control","low_proteinuria","high_proteinuria"))
df$proteinuria_signature=df$immune_impairment_score
stat.test <- df %>%
  t_test(proteinuria_signature ~ proteinuria_cf) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

stat.test <- stat.test %>%
  add_xy_position(x='proteinuria_cf',dodge = 1)

bxp <- ggboxplot(
    df,x="proteinuria_cf",y="proteinuria_signature",
    color= "proteinuria_cf",
    palette=c("pink","blue","green"))

bxp
bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj.signif}",
  tip.length = 0, hide.ns = TRUE
  )

