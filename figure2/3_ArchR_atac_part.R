#!/xtdisk/jiangl_group/huangzh/conda/anaconda3/envs/Enrich1/bin/Rscript

#setwd("F:/Bioinformatics/LiYun/ZhongShanErYuan_breast/UDA_seq/analysis/")

library(ArchR)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pheatmap)
library(data.table)
library(ggplot2)

#register(SnowParam(workers=40,type="SOCK"))
set.seed(1)
addArchRThreads(threads=46)
addArchRGenome("hg38")

donors_batch_1 = fread("/p300s/jiangl_group/ligchA/HuangZheng/kidney/data/intersection_with_scRNA/donors", header=F)
sampleNames_batch_1 = paste0(donors_batch_1$V1, ".batch_1")
donors_batch_2 = fread("/p300s/jiangl_group/ligchA/HuangZheng/kidney/data/intersection_with_scRNA/donors", header=F)
sampleNames_batch_2 = paste0(donors_batch_2$V1, ".batch_2")
sampleNames = c(sampleNames_batch_1, sampleNames_batch_2)

inputFiles_batch_1 = paste0("/p300s/jiangl_group/ligchA/HuangZheng/kidney/data/intersection_with_scRNA/", sampleNames_batch_1, ".atac_fragments.tsv.gz")
inputFiles_batch_2 = paste0("/p300s/jiangl_group/ligchA/HuangZheng/kidney/data/intersection_with_scRNA/", sampleNames_batch_2, ".atac_fragments.tsv.gz")
inputFiles = c(inputFiles_batch_1, inputFiles_batch_2)


# 1.6 Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# 3.1 Creating An ArchRProject
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "projHeme1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# 3.2 Manipulating An ArchRProject
head(projHeme1$cellNames)
head(projHeme1$Sample)
head(projHeme1$TSSEnrichment)
head(projHeme1$ReadsInTSS)
head(projHeme1$PromoterRatio)
head(projHeme1$DoubletScore)
head(projHeme1$DoubletEnrichment)
projHeme1[1:100, ]
projHeme1[projHeme1$cellNames[1:100], ]

#idxPass <- which(projHeme1$TSSEnrichment >= 1.5 & projHeme1$nFrags>=1000)
#projHeme1 <- projHeme1[idxPass, ]

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

# 3.3 Plotting Sample Statistics from an ArchRProject
p1 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)



# 3.4 Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)


# 3.5 Custom plot
cell_count = data.frame(projHeme1$cellNames)
setDT(cell_count)
names(cell_count) = "BARCODE"
cell_count[,donor:=tstrsplit(BARCODE,"\\.",keep=c(1))]
cell_count[,batch_tmp:=tstrsplit(BARCODE,"\\.",keep=c(2))]
cell_count[,batch:=tstrsplit(batch_tmp,"#",keep=c(1))]
cell_count[,group:=ifelse(substr(donor,1,2)=="PB", "disease", "health")]

projHeme1 = addCellColData(
    ArchRProj = projHeme1,
    data  = cell_count$group,
    cells = cell_count$BARCODE,
    name  = "donor_group"
)

meta_info = fread("/p300s/jiangl_group/ligchA/HuangZheng/kidney/analysis/intersection_with_scRNA/scATAC.BARCODE_donor_cellType_UP24h.xls")
meta_info[,donor:=ifelse(substr(donor,1,2)=="WB", substr(donor,1,4), donor)]
#meta_info[,UP24h:=ifelse(UP24h=="Control", 0, UP24h)]
meta_info[,UP24h_group:=ifelse(UP24h=="Control", "UP24h_0g", ifelse(as.double(UP24h)<3, "UP24h_0_to_3g", "UP24h_gt_3g"))]

projHeme1 = addCellColData(
    ArchRProj = projHeme1,
    data  = meta_info$donor,
    cells = meta_info$BARCODE_scATAC,
    name  = "donor"
)

projHeme1 = addCellColData(
    ArchRProj = projHeme1,
    data  = meta_info$predicted.id,
    cells = meta_info$BARCODE_scATAC,
    name  = "cell_type"
)

projHeme1 = addCellColData(
    ArchRProj = projHeme1,
    data  = meta_info$UP24h_group,
    cells = meta_info$BARCODE_scATAC,
    name  = "UP24h_group"
)

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "./projHeme1", load = TRUE)
projHeme1 = loadArchRProject("./projHeme1")

#p1 <- plotFragmentSizes(ArchRProj = projHeme1, groupBy = "donor_group")
p1 <- plotFragmentSizes(ArchRProj = projHeme6, groupBy = "donor_group")
plotPDF(p1, name = "FragSizes.health_VS_disease.pdf", ArchRProj = projHeme6, addDOC = FALSE, width = 5, height = 5)
#p2 <- plotTSSEnrichment(ArchRProj = projHeme1, groupBy = "donor_group")
p2 <- plotTSSEnrichment(ArchRProj = projHeme6, groupBy = "donor_group")
plotPDF(p2, name = "TSS.health_VS_disease.pdf", ArchRProj = projHeme6, addDOC = FALSE, width = 5, height = 5)
#plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.health_VS_disease.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = "QC-Sample-FragSizes-TSSProfile.Sample.pdf", ArchRProj = projHeme6, addDOC = FALSE, width = 5, height = 5)
 
meta = data.frame(projHeme1@cellColData)
setDT(meta, keep.rownames="BARCODE")
fwrite(meta,"meta.intersection.txt",sep="\t")
a = fread("./meta.intersection.txt")
p = ggplot(data=a,aes(x=donor, y=TSSEnrichment)) +
  geom_violin(aes(fill=donor),position="dodge") +
  theme_bw() +
  xlab(element_blank()) +
  ylab("TSS Enrichment") +
  ylim(4,50) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
ggsave("TSS.intersection.png", p, width = 8, height = 3)
ggsave("TSS.intersection.pdf", p, width = 8, height = 3)

p = ggplot(data=a,aes(x=donor, y=log10(nFrags))) +
  geom_violin(aes(fill=donor),position="dodge") +
  theme_bw() +
  xlab(element_blank()) +
  ylab(expression(log[10](Fragments))) +
  ylim(3,5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
ggsave("Fragments.intersection.png", p, width = 8, height = 3)
ggsave("Fragments.intersection.pdf", p, width = 8, height = 3)

# 4.2 Iterative Latent Semantic Indexing (LSI)
projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme1,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 5000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

# 5.1 Clustering using Seurat’s FindClusters() function
projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
head(projHeme2$Clusters)
table(projHeme2$Clusters)
cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM
cM <- cM / Matrix::rowSums(cM)
#p <- pheatmap(
#  mat = as.matrix(cM), 
#  color = paletteContinuous("whiteBlue"), 
#  border_color = "black"
#)
#p

# 6.1 Uniform Manifold Approximation and Projection (UMAP)
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "donor", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "projHeme2", load = TRUE)
projHeme2 = loadArchRProject("./projHeme2")

projHeme3 = loadArchRProject("./projHeme3")
projHeme3 <- addClusters(
  input = projHeme3,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 2,
  maxClusters = 100,
  force = T
)
projHeme3 <- addUMAP(
  ArchRProj = projHeme3, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
p1 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "donor", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.resolution_2.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)


# 9.2 Making Pseudo-bulk Replicates
cell_type_count = data.frame(table(meta_info$predicted.id))
setDT(cell_type_count)
cell_type_count = cell_type_count[Freq>50]
cell_type_gt_50 = as.character(cell_type_count$Var1)
cell_type_gt_50_BARCODE = meta_info[predicted.id %in% cell_type_gt_50, BARCODE_scATAC]
#projHeme3 = subsetCells(ArchRProj = projHeme2, cellNames = cell_type_gt_50_BARCODE)
projHeme3 = subsetArchRProject(ArchRProj=projHeme2, cells=meta_info[predicted.id %in% cell_type_gt_50, BARCODE_scATAC], outputDirectory="projHeme3")
#projHeme3 = loadArchRProject("./projHeme3")

#projHeme3 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters", minReplicates = 2, maxReplicates = 5, sampleRatio = 0.8, minCells=10, maxCells=1000)
projHeme3 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "cell_type", force=T)

# 10.2 peak calling with MACS2 (only linux)
pathToMacs2 <- findMacs2()
projHeme3 <- addReproduciblePeakSet(ArchRProj = projHeme3, groupBy = "cell_type", pathToMacs2 = pathToMacs2, force=T)
getPeakSet(projHeme3)

# 10.4 Add Peak Matrix
projHeme3 <- addPeakMatrix(projHeme3, force=T)
getAvailableMatrices(projHeme3)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
projHeme3 = loadArchRProject("./projHeme3")


p <- plotGroups(
  ArchRProj = projHeme3, 
  groupBy = "donor_group", 
  colorBy = "cellColData", 
  name = "ReadsInPeaks",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + geom_boxplot(width=0.2)
plotPDF(p, name = "ReadsInPeaks.health_vs_disease.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 4, height = 4)

p <- plotGroups(
  ArchRProj = projHeme3, 
  groupBy = "donor", 
  colorBy = "cellColData", 
  name = "ReadsInPeaks",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  ratioYX = 0.3
)
plotPDF(p, name = "ReadsInPeaks.donor.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 10, height = 4)

p <- plotGroups(
  ArchRProj = projHeme3, 
  groupBy = "donor_group", 
  colorBy = "cellColData", 
  name = "FRIP",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p, name = "FRIP.health_vs_disease.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 4, height = 4)

p <- plotGroups(
  ArchRProj = projHeme3, 
  groupBy = "donor", 
  colorBy = "cellColData", 
  name = "FRIP",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  ratioYX = 0.3
)
plotPDF(p, name = "FRIP.donor.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 10, height = 4)

meta = data.frame(projHeme3@cellColData)
setDT(meta, keep.rownames="BARCODE")
fwrite(meta,"./projHeme3/cellColData.txt",sep="\t")
a = fread("./projHeme3/cellColData.txt")
p = ggplot(data=a,aes(x=donor_group, y=ReadsInPeaks)) +
  geom_violin(aes(fill=donor_group),position="dodge", width=0.7) +
  geom_boxplot(outliers = F, width=0.1) +
  theme_bw() +
  xlab(element_blank()) +
  ylab("Reads in peaks") +
  theme(panel.grid = element_blank(),
        legend.position = "none")
ggsave("./projHeme3/Plots/ReadsInPeaks.health_vs_disease.png", p, width = 3, height = 3)
ggsave("./projHeme3/Plots/ReadsInPeaks.health_vs_disease.pdf", p, width = 3, height = 3)

p = ggplot(data=a,aes(x=donor, y=ReadsInPeaks)) +
  geom_violin(aes(fill=donor_group),position="dodge", width=1.3) +
  geom_boxplot(outliers = F, width=0.12) +
  theme_bw() +
  xlab(element_blank()) +
  ylab("Reads in peaks") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
ggsave("./projHeme3/Plots/ReadsInPeaks.donor.png", p, width = 8, height = 3)
ggsave("./projHeme3/Plots/ReadsInPeaks.donor.pdf", p, width = 8, height = 3)

p = ggplot(data=a,aes(x=donor_group, y=FRIP)) +
  geom_violin(aes(fill=donor_group),position="dodge", width=0.7) +
  geom_boxplot(outliers = F, width=0.1) +
  theme_bw() +
  xlab(element_blank()) +
  ylab("FRiP") +
  theme(panel.grid = element_blank(),
        legend.position = "none")
ggsave("./projHeme3/Plots/FRiP.health_vs_disease.png", p, width = 3, height = 3)
ggsave("./projHeme3/Plots/FRiP.health_vs_disease.pdf", p, width = 3, height = 3)

p = ggplot(data=a,aes(x=donor, y=FRIP)) +
  geom_violin(aes(fill=donor_group),position="dodge", width=1) +
  geom_boxplot(outliers = F, width=0.12) +
  theme_bw() +
  xlab(element_blank()) +
  ylab("FRiP") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
p
ggsave("./projHeme3/Plots/FRiP.donor.png", p, width = 8, height = 3)
ggsave("./projHeme3/Plots/FRiP.donor.pdf", p, width = 8, height = 3)

# 11.1 Identifying Marker Peaks with ArchR
#POD = projHeme2[meta_info[predicted.id=="POD", BARCODE_scATAC]]
#markerPeaks <- getMarkerFeatures(
#  ArchRProj = POD, 
#  useMatrix = "PeakMatrix", 
#  groupBy = "UP24h_group",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  testMethod = "wilcoxon"
#)
cell_type_level = fread("./cell_type_level", header=F)
names(cell_type_level) = c("cell_type", "cell_type_level1")
cellColData = fread("./projHeme3/cellColData.txt")
cellColData_1 = merge(cellColData, cell_type_level, by="cell_type")
projHeme3 = addCellColData(
    ArchRProj = projHeme3,
    data  = cellColData_1$cell_type_level1,
    cells = cellColData_1$BARCODE,
    name  = "cell_type_level1"
)
markerPeaks <- getMarkerFeatures(
  ArchRProj = projHeme3, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type_level1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
FDR = SummarizedExperiment::assays(markerPeaks)[["FDR"]]
setDT(FDR)
summary(FDR)
Log2FC = SummarizedExperiment::assays(markerPeaks)[["Log2FC"]]
setDT(Log2FC)
summary(Log2FC)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE) # return a GRanges Object
saveRDS(markerPeaks, "./projHeme3/markerPeaks.celltype_level1.rds")
saveRDS(markerList, "./projHeme3/markerList.celltype_level1.rds")
markerPeaks = readRDS("./projHeme3/markerPeaks.celltype_level1.rds")
markerList = readRDS("./projHeme3/markerList.celltype_level1.rds")

markerGenes <- getMarkerFeatures(
  ArchRProj = projHeme3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "cell_type_level1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerGene_List <- getMarkers(markerGenes, cutOff = "FDR <= 0.01 & Log2FC >= 1") # return a GRanges Object
saveRDS(markerGenes, "./projHeme3/markerGenes.celltype_level1.rds")
saveRDS(markerGene_List, "./projHeme3/markerGene_List.celltype_level1.rds")
markerGenes = readRDS("./projHeme3/markerGenes.celltype_level1.rds")
markerGene_List = readRDS("./projHeme3/markerGene_List.celltype_level1.rds")




# 11.2 Plotting Marker Peaks in ArchR
# 11.2.1 Marker Peak Heatmaps
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = F
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap.celltype_level1", width = 10, height = 30, ArchRProj = projHeme3, addDOC = FALSE)

heatmapGenes_matrix <- plotMarkerHeatmap(
  seMarker = markerGenes, 
  cutOff = "FDR <= 0.1 & Log2FC >= 2",
  transpose = T,
  returnMatrix = T
)
heatmapGenes_matrix = data.frame(t(heatmapGenes_matrix))
setDT(heatmapGenes_matrix, keep.rownames="Gene")


#saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
#projHeme3 = loadArchRProject("./projHeme3")


# 11.2.2 Marker Peak MA and Volcano Plots
pma <- plotMarkers(seMarker = markerPeaks, name = "POD", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pv <- plotMarkers(seMarker = markerPeaks, name = "POD", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "Markers-MA-Volcano.POD.pdf", width = 5, height = 5, ArchRProj = projHeme3, addDOC = FALSE)


# 11.2.3 Marker Peaks in Browser Tracks
#markerGenes  <- c("NEAT1", "CTNND", "ASPH", "NBAS", "FGF13", "NRXN3", "TLL1", "CP", "FOXO3")
p <- plotBrowserTrack(
  ArchRProj = projHeme3, 
  groupBy = "cell_type_level1", 
  region = region,
  #geneSymbol = markerGenes,
  plotSummary = c("bulkTrack", "geneTrack"),
  features = markerList,
  upstream = 50000,
  downstream = 50000,
  tileSize = 1,
  baseSize = 7,
  facetbaseSize = 7
)
plotPDF(plotList = p, name = "Plot-Tracks-With-Features.celltype_level1.pdf", width = 8, height = 5, ArchRProj = projHeme3, addDOC = FALSE)


cell_type_level1 = names(markerGene_List)
for (i in c(1:length(cell_type_level1))) {
    cell_type = cell_type_level1[i]
    cell_type = ifelse(cell_type=="VSM/P", "VSM.P", cell_type)
    cmd_for_eval = paste0("markerGenes_for_track_plot = heatmapGenes_matrix[", cell_type, "==2")
    other_cell_types = c()
    for (j in c(1:length(cell_type_level1))) {
        if (j != i) {
            other_cell_type = cell_type_level1[j]
            other_cell_type = ifelse(other_cell_type=="VSM/P", "VSM.P", other_cell_type)
            cmd_for_eval = paste0(cmd_for_eval, " & ", other_cell_type, "<0.2")
            other_cell_types = c(other_cell_types, other_cell_type)
        }
    }
    cmd_for_eval = paste0(cmd_for_eval, "]")
    eval(parse(text=cmd_for_eval))
    a = markerGenes_for_track_plot[,mget(other_cell_types)]
    a$sd = apply(a,1,sd)
    markerGenes_for_track_plot[,other_cell_types_sd:=a$sd]
    setorder(markerGenes_for_track_plot, other_cell_types_sd)
    if(is.na(markerGenes_for_track_plot[1,Gene])) {print(paste0("No gene passed cutoff in cell: ", cell_type)); next}
    #print(markerGenes_for_track_plot[1:3,Gene])}
    #heatmapGenes.PT = heatmapGenes_matrix[PT==2 & CNT<0 & DCT<0 & DTL<0 & EC<0 & FIB<0 & IC<0 & IMM<0 & PC<0 & PEC<0 & POD<0 & TAL<0 & VSM.P<0][1:3,Gene]
    #markerGenes_for_track_plot = markerGene_List[[i]]$name[1]
    p <- plotBrowserTrack(
        ArchRProj = projHeme3, 
        groupBy = "cell_type_level1", 
        geneSymbol = markerGenes_for_track_plot[1,Gene],
        #geneSymbol = c("CYP4A11", "LCE2D", "MRLN"),
        plotSummary = c("bulkTrack", "geneTrack"),
        features = markerList,
        upstream = 10000,
        downstream = 10000,
        tileSize = 200,
        baseSize = 3,
        facetbaseSize = 7,
        sizes = c(10,2)
    )
    #cell_type = "PT"
    plotPDF(plotList = p, name = paste0("Tracks.markerGenes.celltype_level1.", cell_type,".pdf"), width = 2, height = 5, ArchRProj = projHeme3, addDOC = FALSE)
}



#12.2 Motif Enrichment in Marker Peaks
# devtools::install_github("GreenleafLab/chromVARmotifs") # URL error
projHeme3 <- addMotifAnnotations(ArchRProj = projHeme3, motifSet = "cisbp", name = "Motif", force=T)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
projHeme3 = loadArchRProject("./projHeme3")
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projHeme3,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs
saveRDS(enrichMotifs, "./projHeme3/enrichMotifs.rds")
enrichMotifs = readRDS("./projHeme3/enrichMotifs.rds")


a = fread("./projHeme3/cellColData.txt")
b = data.frame(table(a$cell_type))
setDT(b)
c = b[Freq>600,Var1]
enrichMotifs_main_cell_type = enrichMotifs[,enrichMotifs@colData@rownames %in% c]

heatmapEM <- plotEnrichHeatmap(enrichMotifs_main_cell_type, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "heatmap.Motifs_Enriched_Marker.cell_count_gt_600_per_type", width = 20, height = 14, ArchRProj = projHeme3, addDOC = FALSE)

# 12.3 ArchR Enrichment
# 12.3.1 Encode TF Binding Sites
projHeme3 <- addArchRAnnotations(ArchRProj = projHeme3, collection = "EncodeTFBS", force=T)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
projHeme3 = loadArchRProject("./projHeme3")
enrichEncode <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projHeme3,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichEncode
saveRDS(enrichEncode, "./projHeme3/enrichEncode.rds")
enrichEncode = readRDS("./projHeme3/enrichEncode.rds")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "heatmap.EncodeTFBS_Enriched_Marker", width = 10, height = 14, ArchRProj = projHeme3, addDOC = FALSE)

# 12.3.2 Bulk ATAC-seq
projHeme3 <- addArchRAnnotations(ArchRProj = projHeme3, collection = "ATAC", force = T)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
projHeme3 = loadArchRProject("./projHeme3")
enrichATAC <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projHeme3,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichATAC
saveRDS(enrichATAC, "./projHeme3/enrichATAC.rds")
enrichATAC = readRDS("./projHeme3/enrichATAC.rds")
heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapATAC, name = "heatmap.ATAC_Enriched_Marker", width = 10, height = 14, ArchRProj = projHeme3, addDOC = FALSE)

# 12.2.3 Codex TFBS
projHeme3 <- addArchRAnnotations(ArchRProj = projHeme3, collection = "Codex", force=T)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "projHeme3", load = TRUE)
projHeme3 = loadArchRProject("./projHeme3")
enrichCodex <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projHeme3,
  peakAnnotation = "Codex",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichCodex
saveRDS(enrichCodex, "./projHeme3/enrichCodex.rds")
enrichCodex = readRDS("./projHeme3/enrichCodex.rds")
heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapCodex, name = "heatmap.Codex_Enriched_Marker", width = 10, height = 14, ArchRProj = projHeme3, addDOC = FALSE)

# 13.1 Motif Deviations: POD
POD_BARCODE = meta_info[predicted.id=="POD", BARCODE_scATAC]
#POD = subsetArchRProject(ArchRProj=projHeme3, cells=meta_info[predicted.id %in% POD_BARCODE, BARCODE_scATAC], outputDirectory="POD", force=T)
POD = projHeme3[POD_BARCODE,]
saveArchRProject(ArchRProj = POD, outputDirectory = "POD", load = TRUE)
POD = loadArchRProject("./POD")
POD <- addBgdPeaks(POD)
POD <- addDeviationsMatrix(
  ArchRProj = POD,
  peakAnnotation = "Motif",
  force = TRUE,
  threads = 1
)
saveArchRProject(ArchRProj = POD, outputDirectory = "POD", load = TRUE)
POD = loadArchRProject("./POD")

MotifMatrix = getMatrixFromProject(
  ArchRProj = POD,
  useMatrix = "MotifMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
saveRDS(MotifMatrix, "./POD/MotifMatrix.rds")
MotifMatrix = readRDS("./POD/MotifMatrix.rds")

plotVarDev <- getVarDeviations(POD, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation", width = 5, height = 5, ArchRProj = POD, addDOC = FALSE)



# 13.1 Motif Deviations: EC_GC
EC_GC_BARCODE = meta_info[predicted.id=="EC-GC", BARCODE_scATAC]
#EC_GC = subsetArchRProject(ArchRProj=projHeme3, cells=meta_info[predicted.id %in% EC_GC_BARCODE, BARCODE_scATAC], outputDirectory="EC_GC", force=T)
EC_GC = projHeme3[EC_GC_BARCODE,]
saveArchRProject(ArchRProj = EC_GC, outputDirectory = "EC_GC", load = TRUE)
EC_GC = loadArchRProject("./EC_GC")
EC_GC <- addBgdPeaks(EC_GC)
EC_GC <- addDeviationsMatrix(
  ArchRProj = EC_GC,
  peakAnnotation = "Motif",
  force = TRUE,
  threads = 1
)
saveArchRProject(ArchRProj = EC_GC, outputDirectory = "EC_GC", load = TRUE)
EC_GC = loadArchRProject("./EC_GC")

MotifMatrix = getMatrixFromProject(
  ArchRProj = EC_GC,
  useMatrix = "MotifMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
saveRDS(MotifMatrix, "./EC_GC/MotifMatrix.rds")
MotifMatrix = readRDS("./EC_GC/MotifMatrix.rds")
plotVarDev <- getVarDeviations(EC_GC, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation", width = 5, height = 5, ArchRProj = EC_GC, addDOC = FALSE)
