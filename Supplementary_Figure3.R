
################################################################################
# Supplementary Figure 3: cellchat analysis of cross-species
################################################################################
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)

#input
load("~/Downloads/snRNA/cross_species_snRNA.RData")

Human <- snRNA[,snRNA@meta.data$age %in% c("Human")]
TS <- snRNA[,snRNA@meta.data$age %in% c("TS")]
Macaque <- snRNA[,snRNA@meta.data$age %in% c("Macaque")]
Mouse <- snRNA[,snRNA@meta.data$age %in% c("Mouse")]

{data.input  <- Human@assays[["RNA"]]@data
celltype  <- Human@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = Human$celltype, row.names = names(Human$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = Human$celltype, row.names = names(identify))}

#
{data.input  <- TS@assays[["RNA"]]@data
celltype  <- TS@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = TS$celltype, row.names = names(TS$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = TS$celltype, row.names = names(identify))}

#
{data.input  <- Macaque@assays[["RNA"]]@data
celltype  <- Macaque@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = Macaque$celltype, row.names = names(Macaque$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = Macaque$celltype, row.names = names(identify))
}

#
{data.input  <- Mouse@assays[["RNA"]]@data
  celltype  <- Mouse@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Mouse$celltype, row.names = names(Mouse$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Mouse$celltype, row.names = names(identify))}

##
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat,raw.use = T, population.size = T,type = "truncatedMean",trim = 0.1)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")

ss_Human <- cellchat
saveRDS(ss_Human, "ss_Human.rds")
ss_Macaque <- cellchat
saveRDS(ss_Macaque, "ss_Macaque.rds")
ss_TS <- cellchat
saveRDS(ss_TS, "ss_TS.rds")
ss_Mouse <- cellchat
saveRDS(ss_Mouse, "ss_Mouse.rds")

ss<-data.frame(Human=sort(ss_Human@netP$pathways),Macaque=sort(ss_Macaque@netP$pathways),TS=sort(ss_TS@netP$pathways),Mouse=sort(ss_Mouse@netP$pathways))
Human=data_frame(Human=sort(ss_Human@netP$pathways))
Macaque=data_frame(Macaque=sort(ss_Macaque@netP$pathways))
TS=data_frame(TS=sort(ss_TS@netP$pathways))
Mouse=data_frame(Mouse=sort(ss_Mouse@netP$pathways))
ss <-dplyr::bind_rows(Human,Macaque,TS, Mouse)
write.csv(ss,"cellchat_ss.csv")
ss1<-data.frame(I_A_O_ss=intersect(x=ss_Human@netP$pathways, y = c(ss_Macaque@netP$pathways, ss_TS@netP$pathways,ss_Mouse@netP$pathways)))
write.csv(ss1,"cellchat_ss1.csv")

ECM_Human <- cellchat
saveRDS(ECM_Human, "ECM_Human.rds")
ECM_Macaque <- cellchat
saveRDS(ECM_Macaque, "ECM_Macaque.rds")
ECM_TS <- cellchat
saveRDS(ECM_TS, "ECM_TS.rds")
ECM_Mouse <- cellchat
saveRDS(ECM_Mouse, "ECM_Mouse.rds")
ECM<-data.frame(Human=sort(ECM_Human@netP$pathways),Macaque=sort(ECM_Macaque@netP$pathways),TS=sort(ECM_TS@netP$pathways),Mouse=sort(ECM_Mouse@netP$pathways))
write.csv(ECM,"cellchat_ECM.csv") 
Human=data_frame(Human=sort(ECM_Human@netP$pathways))
Macaque=data_frame(Macaque=sort(ECM_Macaque@netP$pathways))
TS=data_frame(TS=sort(ECM_TS@netP$pathways))
Mouse=data_frame(Mouse=sort(ECM_Mouse@netP$pathways))

ECM <-dplyr::bind_rows(Human,Macaque,TS, Mouse)
write.csv(ECM,"cellchat_ECM.csv")
ECM1<-data.frame(I_A_O_ECM=intersect(x=ECM_Human@netP$pathways, y = c(ECM_Macaque@netP$pathways, ECM_TS@netP$pathways, ECM_Mouse@netP$pathways)))
write.csv(ECM1,"cellchat_ECM1.csv")

ccc_Human <- cellchat
saveRDS(ccc_Human, "ccc_Human.rds")
ccc_Macaque <- cellchat
saveRDS(ccc_Macaque, "ccc_Macaque.rds")
ccc_TS <- cellchat
saveRDS(ccc_TS, "ccc_TS.rds")
ccc_Mouse <- cellchat
saveRDS(ccc_Mouse, "ccc_Mouse.rds")
ccc<-data.frame(Human=sort(ccc_Human@netP$pathways),Macaque=sort(ccc_Macaque@netP$pathways),TS=sort(ccc_TS@netP$pathways),Mouse=sort(ccc_Mouse@netP$pathways))
write.csv(ccc,"cellchat_ccc.csv") 
ccc<-data.frame(Human=sort(ccc_Human@netP$pathways),Macaque=sort(ccc_Macaque@netP$pathways),TS=sort(ccc_TS@netP$pathways),Mouse=sort(ccc_Mouse@netP$pathways))
Human=data_frame(Human=sort(ccc_Human@netP$pathways))
Macaque=data_frame(Macaque=sort(ccc_Macaque@netP$pathways))
TS=data_frame(TS=sort(ccc_TS@netP$pathways))
Mouse=data_frame(Mouse=sort(ccc_Mouse@netP$pathways))
ccc <-dplyr::bind_rows(Human,Macaque,TS, Mouse)
write.csv(ccc,"cellchat_ccc.csv") 
ccc1<-data.frame(I_A_O_ccc=intersect(x=ccc_Human@netP$pathways, y = c(ccc_Macaque@netP$pathways, ccc_TS@netP$pathways,ccc_Mouse@netP$pathways)))
write.csv(ccc1,"cellchat_ccc1.csv") 

# view results
cellchat@netP$pathways
head(cellchat@LR$LRsig)

# merge cellchat object
con.list <- list(Human=ss_Human,Macaque=ss_Macaque,TS=ss_TS,Mouse=ss_Mouse)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
saveRDS(cellchat, "cellchat_ss.rds")
con.list <- list(Human=ECM_Human,Macaque=ECM_Macaque,TS=ECM_TS,Mouse=ECM_Mouse)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
saveRDS(cellchat, "cellchat_ECM.rds")
con.list <- list(Human=ccc_Human,Macaque=ccc_Macaque,TS=ccc_TS,Mouse=ccc_Mouse)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
saveRDS(cellchat, "cellchat_ccc.rds")

################################################################################
# Supplementary Figure 3A
################################################################################
con.list[[1]]@net$count
con.list[[2]]@net$count
con.list[[3]]@net$count
con.list[[4]]@net$count
# Prism 9 visualizes the number of cell communications for each cell type of each species

################################################################################
# Supplementary Figure 3B
# Identification and visualization of conserved and specific signaling pathways
################################################################################
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2, 3, 4),stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2, 3, 4),stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave("compare_pathway_strength.pdf", p,width = 10,height = 8)

################################################################################
# Supplementary Figure 3C-E
# Representative pathways with clear differences between species
################################################################################

{pdf("Compare_signal_net.pdf", width = 12,height = 5.5)
pathways.show <-c("EGF")
weight.max <- getMaxWeight(con.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,4), xpd=T)
for (i in 1:length(con.list)) {
  netVisual_aggregate(con.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 8,
                      vertex.size.max = 7,vertex.size = groupSize,
                      signaling.name = paste(pathways.show, names(con.list)[i]))  
}

pathways.show <-c("PTN")
weight.max <- getMaxWeight(con.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,4), xpd=T)
for (i in 1:length(con.list)) {
  netVisual_aggregate(con.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 8,
                      vertex.size.max = 7,vertex.size = groupSize,
                      signaling.name = paste(pathways.show, names(con.list)[i]))  
}
pathways.show <-c("IGF")
weight.max <- getMaxWeight(con.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,4), xpd=T)
for (i in 1:length(con.list)) {
  netVisual_aggregate(con.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 8,
                      vertex.size.max = 7,vertex.size = groupSize,
                      signaling.name = paste(pathways.show, names(con.list)[i]))  
}

#
{pdf("Compare_signal_net_2.pdf", width = 12,height = 5.5)
  pathways.show <-c("VISTA")
  weight.max <- getMaxWeight(con.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,4), xpd=T)
  for (i in 1:length(con.list)) {
    netVisual_aggregate(con.list[[i]],signaling = pathways.show, layout = "circle",
                        edge.weight.max = weight.max[1], edge.width.max = 8,
                        vertex.size.max = 7,vertex.size = groupSize,
                        signaling.name = paste(pathways.show, names(con.list)[i]))  
  }
  dev.off()
}

