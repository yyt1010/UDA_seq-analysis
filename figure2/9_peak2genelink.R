library(ArchR)
library(parallel)
addArchRThreads(threads = 6, force = FALSE)

data=readRDS("add_gene.rds")  ##added gene expression matrix in ArchR object
data <- addCoAccessibility(
    ArchRProj = data,
    reducedDims = "IterativeLSI"
)

data <- addPeak2GeneLinks(
    ArchRProj = data,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix"
)
p1 <- plotEmbedding(ArchRProj = data, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
p1
p1 <- plotPeak2GeneHeatmap(ArchRProj = data, groupBy = "celltype")
print(p1)
plotPDF(p1,
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks3.pdf",
    ArchRProj = data,
    addDOC = FALSE)

saveRDS(data,file="9_addCo.rds")

