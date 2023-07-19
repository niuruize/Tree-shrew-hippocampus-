
################################################################################
# Supplementary Figure 7
# Analysis of cell communication during TS hippocampal aging
################################################################################
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)

#input data
Infancy<- AB[,AB@meta.data$age %in% c("Infancy")]
Adult <- AB[,AB@meta.data$age %in% c("Adult")]
Old <- AB[,AB@meta.data$age %in% c("Old")]

{data.input  <- Infancy@assays[["RNA"]]@data
celltype  <- Infancy@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = Infancy$celltype, row.names = names(Infancy$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = Infancy$celltype, row.names = names(identify))}
##
{data.input  <- Old@assays[["RNA"]]@data
celltype  <- Old@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = Old$celltype, row.names = names(Old$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = Old$celltype, row.names = names(identify))}
##
{data.input  <- Adult@assays[["RNA"]]@data
celltype  <- Adult@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = Adult$celltype, row.names = names(Adult$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = Adult$celltype, row.names = names(identify))
}

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
# set the used database in the object
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
######plan("multiprocess", workers = 4) # do parallel  
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

ss_Infancy <- cellchat
saveRDS(ss_Infancy, "ss_Infancy.rds")
ss_Adult <- cellchat
saveRDS(ss_Adult, "ss_Adult.rds")
ss_Old <- cellchat
saveRDS(ss_Old, "ss_Old.rds")

ss<-data.frame(Infancy=sort(ss_Infancy@netP$pathways),Adult=sort(ss_Adult@netP$pathways),Old=sort(ss_Old@netP$pathways))
Infancy=data_frame(Infancy=sort(ss_Infancy@netP$pathways))
Adult=data_frame(Adult=sort(ss_Adult@netP$pathways))
Old=data_frame(Old=sort(ss_Old@netP$pathways))
ss <-dplyr::bind_rows(Infancy,Adult,Old)
write.csv(ss,"cellchat_ss.csv") 
ss1<-data.frame(I_A_O_ss=intersect(x=ss_Infancy@netP$pathways, y = c(ss_Adult@netP$pathways, ss_Old@netP$pathways)))
write.csv(ss1,"cellchat_ss1.csv") 

ECM_Infancy <- cellchat
saveRDS(ECM_Infancy, "ECM_Infancy.rds")
ECM_Adult <- cellchat
saveRDS(ECM_Adult, "ECM_Adult.rds")
ECM_Old <- cellchat
saveRDS(ECM_Old, "ECM_Old.rds")
ECM<-data.frame(Infancy=sort(ECM_Infancy@netP$pathways),Adult=sort(ECM_Adult@netP$pathways),Old=sort(ECM_Old@netP$pathways))
write.csv(ECM,"cellchat_ECM.csv") 
Infancy=data_frame(Infancy=sort(ECM_Infancy@netP$pathways))
Adult=data_frame(Adult=sort(ECM_Adult@netP$pathways))
Old=data_frame(Old=sort(ECM_Old@netP$pathways))
ECM <-dplyr::bind_rows(Infancy,Adult,Old)
write.csv(ECM,"cellchat_ECM.csv") 
ECM1<-data.frame(I_A_O_ECM=intersect(x=ECM_Infancy@netP$pathways, y = c(ECM_Adult@netP$pathways, ECM_Old@netP$pathways)))
write.csv(ECM1,"cellchat_ECM1.csv")

ccc_Infancy <- cellchat
saveRDS(ccc_Infancy, "ccc_Infancy.rds")
ccc_Adult <- cellchat
saveRDS(ccc_Adult, "ccc_Adult.rds")
ccc_Old <- cellchat
saveRDS(ccc_Old, "ccc_Old.rds")
ccc<-data.frame(Infancy=sort(ccc_Infancy@netP$pathways),Adult=sort(ccc_Adult@netP$pathways),Old=sort(ccc_Old@netP$pathways))
write.csv(ccc,"cellchat_ccc.csv")

ccc<-data.frame(Infancy=sort(ccc_Infancy@netP$pathways),Adult=sort(ccc_Adult@netP$pathways),Old=sort(ccc_Old@netP$pathways))
Infancy=data_frame(Infancy=sort(ccc_Infancy@netP$pathways))
Adult=data_frame(Adult=sort(ccc_Adult@netP$pathways))
Old=data_frame(Old=sort(ccc_Old@netP$pathways))
ccc <-dplyr::bind_rows(Infancy,Adult,Old)
write.csv(ccc,"cellchat_ccc.csv") 
ccc1<-data.frame(I_A_O_ccc=intersect(x=ccc_Infancy@netP$pathways, y = c(ccc_Adult@netP$pathways, ccc_Old@netP$pathways)))
write.csv(ccc1,"cellchat_ccc1.csv")

# 
con.list <- list(Infancy=ss_Infancy,Adult=ss_Adult,Old=ss_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
con.list <- list(Infancy=ECM_Infancy,Adult=ECM_Adult,Old=ECM_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
con.list <- list(Infancy=ccc_Infancy,Adult=ccc_Adult,Old=ccc_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)

################################################################################
# Supplementary Figure 7 A
# Identification and visualization of conserved and specific signaling pathways
################################################################################
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2, 3),stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1, 2, 3),stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave("compare_pathway_strength.pdf", p,width = 10,height = 5.5)

################################################################################
# Supplementary Figure 7B
# Overall view of all cell populations: number versus intensity of communication
################################################################################
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = c("#91B822","#FFE76F","#CDA59E") ,group = c(1,2,3),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, color.use = c("#91B822","#FFE76F","#CDA59E") ,group = c(1,2,3),measure = "weight")
p <- gg1+gg2
ggsave("overview_number_strength.pdf", p, width = 6,height = 4)
dev.off()

################################################################################
# Supplementary Figure 7C
# Overall view of all cell populations: number versus intensity of communication
################################################################################
con.list[[1]]@net$coun
con.list[[2]]@net$count
con.list[[3]]@net$count

pdf("Compare_signal_net.pdf", width = 10,height = 10)
par(mfrow = c(2,2))
h1 <- netVisual_diffInteraction(cellchat,edge.weight.max = 10,comparison = c(1, 2), weight.scale = T)
h2 <- netVisual_diffInteraction(cellchat, measure = "weight",comparison = c(1, 2),vertex.size.max = 15, weight.scale = T)
h3 <- netVisual_diffInteraction(cellchat,edge.weight.max = 10,comparison = c(2, 3), weight.scale = T)
h4 <- netVisual_diffInteraction(cellchat, measure = "weight",comparison = c(2, 3),vertex.size.max = 15, weight.scale = T)
dev.off()

