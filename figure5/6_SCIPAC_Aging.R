setwd("/p300s/jiangl_group/chenyanjie/research/UDA/Step6_SCIPAC")
library(SCIPAC)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggthemes)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

sc.dat <- readRDS("/p300s/jiangl_group/huangzh/Post_index/PBMC/UDA/QC_clusters/TCR/l2_new.rds")
counts_matrix <- sc.dat@assays$RNA@counts
ordinary_matrix <- as.matrix(counts_matrix)
print(ordinary_matrix[1:100, 1:10])
sc.meta <- sc.dat@meta.data

# Check the first few rows
head(sc.meta)
table(sc.meta$donor)
DefaultAssay(sc.dat) <- "RNA"
human_sce <- as.SingleCellExperiment(sc.dat)
pb <- aggregateData(human_sce,
                    assay = "counts", 
                    fun = "sum",
                    by = c("donor"))
bulkmatrix <- pb@assays@data[[1]]
dge <- DGEList(counts = bulkmatrix)
dge <- calcNormFactors(dge)
normalized_counts_tmm <- cpm(dge, normalized.lib.sizes = TRUE)
saveRDS(normalized_counts_tmm,"normalized_counts_tmm.RDS")

bulk.dat <- readRDS("normalized_counts_tmm.RDS")
print(bulk.dat[1:4, 1:6])
bulk.meta <- unique(sc.meta[, c("donor", "age")])
rownames(bulk.meta)<-bulk.meta[,1]
colnames(bulk.meta)<-c("case.id","phenotype")
head(bulk.meta)

obtain.preprocessed.data_modi <- function(exprs.data, hvg = 1000){
  exprs.data <- Seurat::CreateSeuratObject(exprs.data, project = "CreateSeuratObject", min.cells = 3, min.features = 200)
  exprs.data <- Seurat::NormalizeData(exprs.data, normalization.method = "LogNormalize", scale.factor = 10000)
  exprs.data <- Seurat::FindVariableFeatures(exprs.data, selection.method = "vst", nfeatures = hvg)
  exprs.data[['RNA']] = as(object = exprs.data[['RNA']], Class = "Assay") # Add
  exprs.data.norm <- exprs.data[["RNA"]]@data
  exprs.data.variable.genes <- Seurat::VariableFeatures(exprs.data)
  processed.data <- exprs.data.norm[exprs.data.variable.genes, ] %>%
    as.matrix()
  return(processed.data)
}

preprocess.sc.bulk.dat_modi <- function(sc.dat, bulk.dat, hvg = 1000){
  overlap.genes <- intersect(rownames(sc.dat), rownames(bulk.dat))
  sc.dat.new <- sc.dat[overlap.genes, ]
  bulk.dat.new <- bulk.dat[overlap.genes, ] %>% as.matrix()
  bulk.dat.new <- log(bulk.dat.new + 1)
  sc.dat.preprocessed <- obtain.preprocessed.data_modi(exprs.data = sc.dat.new, hvg = hvg)
  bulk.dat.preprocessed <- bulk.dat.new[rownames(sc.dat.preprocessed), ]
  return(list("sc.dat.preprocessed" = sc.dat.preprocessed,
              "bulk.dat.preprocessed" = bulk.dat.preprocessed))
}

sc.bulk.prec <- preprocess.sc.bulk.dat_modi(ordinary_matrix, bulk.dat, hvg = 1000)
sc.dat.prep <- sc.bulk.prec$sc.dat.preprocessed
bulk.dat.prep <- sc.bulk.prec$bulk.dat.preprocessed
dim(sc.dat.prep)
dim(bulk.dat.prep)
pca.res <- sc.bulk.pca(sc.dat.prep, bulk.dat.prep, do.pca.sc = FALSE, n.pc = 30)
sc.dat.dim.redu <- pca.res$sc.dat.rot
bulk.dat.dim.redu <- pca.res$bulk.dat.rot
ct.res <- seurat.ct(sc.dat.dim.redu, res = 1.5)
summary(ct.res)

k <- ct.res$k
print(k)
# [1] 28

ct.assignment <- ct.res$ct.assignment
head(ct.assignment)
centers <- ct.res$centers
print(centers[1:6, 1:3])

classifier.Lambda	 <- function(bulk.dat, y, family, K.means.res, ela.net.alpha, 
                              bt.size, nfold = nfold, numCores){
  fx <- function(seed){
    set.seed(seed)
    return(classifier.Lambda.core(bulk.dat, y, family, K.means.res, 
                                  ela.net.alpha = ela.net.alpha, nfold = nfold))
  }
  K <- K.means.res$k
  ct.assign <- K.means.res$ct.assignment

  seed.ls <- c(1:bt.size)
  Lambda.tab <- parallel::mclapply(seed.ls, fx, mc.cores = numCores)

  Lambda.res <- matrix(NA, nrow = nrow(ct.assign), ncol = bt.size)
  for (i in 1:ncol(Lambda.res)) {
    Lambda.res[, i] <- Lambda.tab[[i]]$Lambda
  }

  # Delete columns with NAs.
  if(sum(is.na(Lambda.res[1, ])) == 0){
    Lambda.res <- Lambda.res
  } else {
    na.idx <- which(is.na(Lambda.res[1, ]))
    Lambda.res <- Lambda.res[, -na.idx]
  }
  return(Lambda.res)
}

bulk.meta.new <- bulk.meta[rownames(bulk.dat.dim.redu),]
y <- bulk.meta.new[,2]

SCIPAC.res <- SCIPAC(bulk.dat = bulk.dat.dim.redu, y = y, 
                     family = "gaussian", ct.res = ct.res, 
                     ela.net.alpha = 0.4, bt.size = 50, 
                     numCores = 40, CI.alpha = 0.05, nfold = 10)

head(SCIPAC.res)
saveRDS(SCIPAC.res,"SCIPAC.res.RDS")

LUAD.sc.umap <- sc.dat@reductions$umap@cell.embeddings

LUAD.sc.meta <- sc.meta[rownames(SCIPAC.res), ]
all(rownames(sc.meta) == rownames(SCIPAC.res))
plot.dat <- cbind(sc.meta, SCIPAC.res)
all(rownames(plot.dat) == rownames(LUAD.sc.umap))
plot.dat <- cbind(plot.dat, LUAD.sc.umap)
plot.dat$log.pval.adj <- ifelse((plot.dat$sig == "Not.sig"),
                                0, plot.dat$log.pval)

# Plot the cell types
colors.cell <- tableau_color_pal("Classic Cyclic",
                                 direction = 1)(length(unique(plot.dat$predicted.celltype.l2)))
plt <- plot.dat %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, col = predicted.celltype.l2, fill = predicted.celltype.l2)) +
  theme_tufte(base_size = 12) +
  theme(panel.background = element_rect(fill = NA, color = "black")) +
  guides(color = guide_legend(override.aes = list(stroke = 1,
                                                  alpha = 1, shape = 16, size = 4)),
         alpha = "none") +
  scale_color_manual(" ", values = colors.cell, guide = "none") +
  scale_fill_manual(" ", values = colors.cell, guide = "none") +
  theme(legend.position= "top") +
  labs(x = "UMAP 1", y = "UMAP 2")


p1 <- plt + geom_point()

# Control the range of the estimate of Lambda.Psi within -2 and 2
plot.dat$Lambda.adj <- plot.dat$Lambda.est
if(sum(plot.dat$Lambda.est <= -2) > 0){
  plot.dat[which(plot.dat$Lambda.est <= -2), ]$Lambda.adj <- -2
}
if(sum(plot.dat$Lambda.est >= 2) > 0){
  plot.dat[which(plot.dat$Lambda.est >= 2), ]$Lambda.adj <- 2
}

p2 <- ggplot(plot.dat) + 
  geom_point(aes(x=UMAP_1, y=UMAP_2, col=Lambda.adj),, size=0.2) +
  scale_color_gradient2(expression(Lambda), low="darkblue",high='red', mid = "white", midpoint = 0,
                        limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2),
                        labels = c(expression(phantom(x) <= -2), -1, 0, 1, expression(phantom(x) >= 2))) +
  theme(legend.title = element_text(family = "sans", size = 27),
        legend.position= "top",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "lightgrey", linewidth = 0.1))+ 
  labs(x = "UMAP 1", y = "UMAP 2")



plot.dat$pval.adj.adj <- plot.dat$log.pval.adj
if(sum(plot.dat$log.pval.adj <= -3) > 0){
  plot.dat[which(plot.dat$log.pval.adj <= -3), ]$pval.adj.adj <- -3
}
if(sum(plot.dat$log.pval.adj >= 3) > 0){
  plot.dat[which(plot.dat$log.pval.adj >= 3), ]$pval.adj.adj <- 3
}

plot.dat.positive <- plot.dat[which(plot.dat$pval.adj.adj > 0), ]
plot.dat.positive$pval.adj.adj <- -plot.dat.positive$pval.adj.adj

plot.dat.negative <- plot.dat[which(plot.dat$log.pval.adj <= 0), ]

plot.dat.negative <- plot.dat.negative %>% arrange(desc(pval.adj.adj))
plot.dat.positive <- plot.dat.positive %>% arrange(desc(pval.adj.adj))

# Plot
p3 <- ggplot() + 
  geom_point(data = plot.dat.negative, aes(x=UMAP_1, y=UMAP_2, col=pval.adj.adj), size=0.2) +
  scale_color_gradient2('   ', low="darkblue", high='white', mid="white", midpoint=0,
                        limits=c(-3, -1.3), breaks=c(-3, -2, -1.3), 
                        labels=c(expression(phantom(x) <= 0 .001), 0.01, 0.05),
                        na.value="lightgrey") +
  new_scale_color() + 
  geom_point(data = plot.dat.positive, aes(x=UMAP_1, y=UMAP_2, col=pval.adj.adj), size=0.2) +
  scale_color_gradient2("p-values", low="red", high='white', mid="white", midpoint=0,
                        limits=c(-3, -1.3), breaks=c(-3, -2, -1.3), 
                        labels=c(expression(phantom(x) <= 0.001), 0.01, 0.05),
                        na.value="lightgrey") + 
  theme(legend.position="top",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_line(color="lightgrey", size=0.1)) + 
  labs(x="UMAP 1", y="UMAP 2")

pdf("SCIPAC.pdf",width=20,height=5)
cowplot::plot_grid(p1, p2, p3, nrow = 1)
dev.off()

pdf("SCIPAC_2.pdf",width=7,height=7)
p2
dev.off()

pdf("SCIPAC_3.pdf",width=7,height=7)
p3
dev.off()


# seurat_clusters
Idents(sc.dat)<-sc.dat$SCIPAC_label
markers_age_pos_none <- FindMarkers(sc.dat, ident.1 = "Age_pos", ident.2 = "None")

markers_age_neg_none <- FindMarkers(sc.dat, ident.1 = "Age_neg", ident.2 = "None")


dif_pos=data.frame(
  symbol=rownames(markers_age_pos_none),
  log2FoldChange=markers_age_pos_none$avg_log2FC,
  padj=markers_age_pos_none$p_val_adj
)
# dif_pos <- dif_pos[dif_pos$log2FoldChange >= 1, ]
# dif_pos <- dif_pos %>% arrange(dif_pos$padj)

dim(dif_pos)
# [1] 424   3
write.csv(dif_pos,"Age_pos_Vs_None.csv")

dif_neg=data.frame(
  symbol=rownames(markers_age_neg_none),
  log2FoldChange=markers_age_neg_none$avg_log2FC,
  padj=markers_age_neg_none$p_val_adj
)
# dif_neg <- dif_neg[dif_neg$log2FoldChange >=1,]
# dif_neg <- dif_neg %>% arrange(dif_neg$padj)

dim(dif_neg)
# [1] 339   3
write.csv(dif_neg,"Age_neg_Vs_None.csv")


colors <- c("#E07B91", "#F0B98D", "#C7C7C7", "black")

age_label <- sc.dat@meta.data$SCIPAC_label
cell_type <- sc.dat@meta.data$predicted.celltype.l2

data <- data.frame(age_label, cell_type)

cell_counts_pos <- table(data$cell_type[data$age_label %in% c("Age_pos")])

plot_data <- data.frame(
  cell_type = names(cell_counts_pos),
  cell_count = as.vector(cell_counts_pos)
)

plot_data$percentage <- plot_data$cell_count / sum(plot_data$cell_count) * 100

plot_data <- plot_data %>% arrange(desc(percentage))

plot_data$cell_type[4:nrow(plot_data)] <- "Other"

plot_data <- plot_data %>%
  group_by(cell_type) %>%
  summarise(
    cell_count = sum(cell_count),
    percentage = sum(percentage)
  ) %>%
  ungroup()

plot_data$legend_label <- paste0(plot_data$cell_type, " (", round(plot_data$percentage, 2), "%)")

pie_plot <- ggplot(plot_data, aes(x = "", y = cell_count, fill = factor(legend_label, levels = plot_data$legend_label))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  labs(
    fill = "",
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_manual(values = colors) +
  geom_text(data = plot_data %>% filter(cell_type != "Other"), aes(label = cell_type), position = position_stack(vjust = 0.5), size = 4)

pdf("pie_plot_pos.pdf", width = 6, height = 6)
print(pie_plot)
dev.off()

colors <- c("#8595E1", "#D6BCC0", "#8E063B", "black")

age_label <- sc.dat@meta.data$SCIPAC_label
cell_type <- sc.dat@meta.data$predicted.celltype.l2

data <- data.frame(age_label, cell_type)

cell_counts_neg <- table(data$cell_type[data$age_label %in% c("Age_neg")])

plot_data <- data.frame(
  cell_type = names(cell_counts_neg),
  cell_count = as.vector(cell_counts_neg)
)

plot_data$percentage <- plot_data$cell_count / sum(plot_data$cell_count) * 100

plot_data <- plot_data %>% arrange(desc(percentage))

plot_data$cell_type[4:nrow(plot_data)] <- "Other"

plot_data <- plot_data %>%
  group_by(cell_type) %>%
  summarise(
    cell_count = sum(cell_count),
    percentage = sum(percentage)
  ) %>%
  ungroup()

plot_data$legend_label <- paste0(plot_data$cell_type, " (", round(plot_data$percentage, 2), "%)")

pie_plot <- ggplot(plot_data, aes(x = "", y = cell_count, fill = factor(legend_label, levels = plot_data$legend_label))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  labs(
    fill = "",
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_manual(values = colors) +
  geom_text(data = plot_data %>% filter(cell_type != "Other"), aes(label = cell_type), position = position_stack(vjust = 0.5), size = 4, hjust = 0.5)

pdf("pie_plot_neg.pdf", width = 6, height = 6)
print(pie_plot)
dev.off()
