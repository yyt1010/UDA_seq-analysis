setwd("/p300s/jiangl_group/chenyanjie/research/UDA/Step5_CCC/02_UP24h/output")

library(CellChat)

########################################
# compare multiple datasets
########################################
# AAV
data.input <- readRDS("/xtdisk/jiangl_group/chenyanjie/research/UDA/Step5_CCC/03_POD_EC-GC/output/1_mat_AAV.rds")
meta <- readRDS("/xtdisk/jiangl_group/chenyanjie/research/UDA/Step5_CCC/03_POD_EC-GC/output/1_meta_AAV.rds")
# [1] 36601  4649


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

saveRDS(cellchat,file="cellchat_AAV.rds")

# Control
data.input <- readRDS("/xtdisk/jiangl_group/chenyanjie/research/UDA/Step5_CCC/03_POD_EC-GC/output/1_mat_Control.rds")
meta <- readRDS("/xtdisk/jiangl_group/chenyanjie/research/UDA/Step5_CCC/03_POD_EC-GC/output/1_meta_Control.rds")
# [1] 36601 3643

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

saveRDS(cellchat,file="cellchat_Control.rds")

# AAV_vs_Control
setwd("3g_vs_Control/")

cellchat.AAV <- readRDS("/p300s/jiangl_group/chenyanjie/research/UDA/Step5_CCC/02_UP24h/output/cellchat_AAV.rds")
cellchat.Control <- readRDS("/p300s/jiangl_group/chenyanjie/research/UDA/Step5_CCC/02_UP24h/output/cellchat_Control.rds")


object.list <- list(AAV = cellchat.AAV, Control = cellchat.Control)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#saveRDS(cellchat,file="cellchat_merge.rds")

pdf("compare1.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

cellchat.AAV@netP$pathways

cellchat.Control@netP$pathways

object.list <- list(AAV = cellchat.AAV, Control = cellchat.Control)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf("compare1.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

 # 744 signaling genes.
 # 49196 cells. 

pdf("compare2.pdf",width=5, height=5)
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


pdf("compare3.pdf", width=10, height=5)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()


pdf("compare11.pdf",width=10, height=5)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()


