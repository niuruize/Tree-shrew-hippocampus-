
################################################################################
# Figure 8
# Endo
################################################################################
library(Seurat)
library(tidyverse)
library(liger)
library(patchwork)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(Nebulosa) 

################################################################################
# Figure 8A
################################################################################
load("~/Downloads/snRNA/snRNA.RData")
#Wnt_score
Wnt_gene <- read_xlsx("Wnt.xlsx")
gene <- as.list(Wnt_gene)
snRNA <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Wnt")
colnames(snRNA@meta.data)[11]<-"Wnt_Score" 
##density
tmp= snRNA[,snRNA@meta.data$celltype %in% c("Endo")]
b<-FetchData(tmp, vars = c("orig.ident","Wnt_Score","celltype"))
b[["orig.ident"]]<-factor(b[["orig.ident"]], levels=c("Infancy","Adult","Old"))
P1<-ggdensity(b, x = "Wnt_Score",
              add = "mean", rug = TRUE,
              color = "orig.ident", fill = "orig.ident",
              palette = c("#00AFBB", "#E7B800","#8B658B"))
ggsave(filename = "Wnt_Score.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 

################################################################################
# Figure 8B, C
################################################################################
load("~/Downloads/snRNA/snRNA.RData")
snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("Endo")]
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 10, verbose = FALSE)
snRNA <- FindNeighbors(snRNA, dims = 1:10)
snRNA <- FindClusters(snRNA, resolution = 0.1)
snRNA <- RunUMAP(snRNA, dims = 1:10)
snRNA <- RunTSNE(snRNA, dims = 1:10)
head(Idents(snRNA), 5) 
current.cluster.ids <- c("0","1","2")
new.cluster.ids <- c("End1","End2","End3") 
snRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
snRNA <- RenameIdents(snRNA, `0`="End1",`1`="End2",`2`="End3")
Idents(snRNA) <- factor(Idents(snRNA), levels=c("End1","End2","End3"))
count_table <- table(snRNA@meta.data[["seurat_clusters"]], snRNA@meta.data[["age"]])
count_table
save(snRNA,file="TS_Endo.RData")
p1 <- DimPlot(snRNA, reduction = "umap", group.by = "age", pt.size=1)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p2 <- DimPlot(snRNA, reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap1 <- plot_grid(p1, p2,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 24, height = 9, units = 'cm')

################################################################################
# Figure 8D
################################################################################
Wnt_gene<-read_xlsx("Wnt.xlsx")
gene<-as.list(Wnt_gene)
snRNA<-AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Wnt")
colnames(snRNA@meta.data)[13]<-"Wnt_Score"
#boxplot
b<-FetchData(snRNA, vars = c("orig.ident","Wnt_Score","celltype"))
ggboxplot(b, x = "orig.ident", y = "Wnt_Score",facet.by = "celltype",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "orig.ident", palette = "npg",ncol=3)+stat_compare_means(comparisons = my_comparisons)
ggsave("Endo_subcelltype.pdf",width = 20,height = 10,units = "cm")

################################################################################
# Figure 8H,I
# cellchat
################################################################################
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
#input
Infancy<- snRNA[,snRNA@meta.data$age %in% c("Infancy")]
Adult <- snRNA[,snRNA@meta.data$age %in% c("Adult")]
Old <- snRNA[,snRNA@meta.data$age %in% c("Old")]

{data.input  <- Infancy@assays[["RNA"]]@data
  celltype  <- Infancy@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Infancy$celltype, row.names = names(Infancy$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Infancy$celltype, row.names = names(identify))}

{data.input  <- Old@assays[["RNA"]]@data
  celltype  <- Old@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Old$celltype, row.names = names(Old$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Old$celltype, row.names = names(identify))}

{data.input  <- Adult@assays[["RNA"]]@data
  celltype  <- Adult@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Adult$celltype, row.names = names(Adult$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Adult$celltype, row.names = names(identify))
}

#
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

#merge
con.list <- list(Infancy=ss_Infancy,Adult=ss_Adult,Old=ss_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
con.list <- list(Infancy=ECM_Infancy,Adult=ECM_Adult,Old=ECM_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
con.list <- list(Infancy=ccc_Infancy,Adult=ccc_Adult,Old=ccc_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
saveRDS(cellchat, "cellchat.rds")

#
par(mfrow = c(3,1))
s.cell <- c("Microglia","End","MOL")

con.list <- list(Infancy=ss_Infancy,Adult=ss_Adult,Old=ss_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
count_ss1 <- con.list[[1]]@net$count[s.cell, s.cell]
count_ss2 <- con.list[[2]]@net$count[s.cell, s.cell]
count_ss3 <- con.list[[3]]@net$count[s.cell, s.cell]

con.list <- list(Infancy=ECM_Infancy,Adult=ECM_Adult,Old=ECM_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
count_ECM1 <- con.list[[1]]@net$count[s.cell, s.cell]
count_ECM2 <- con.list[[2]]@net$count[s.cell, s.cell]
count_ECM3 <- con.list[[3]]@net$count[s.cell, s.cell]

con.list <- list(Infancy=ccc_Infancy,Adult=ccc_Adult,Old=ccc_Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
count_cc1 <- con.list[[1]]@net$count[s.cell, s.cell]
count_cc2 <- con.list[[2]]@net$count[s.cell, s.cell]
count_cc3 <- con.list[[3]]@net$count[s.cell, s.cell]

count1 <- count_ss1+count_cc1+count_ECM1
count2 <- count_ss2+count_cc2+count_ECM2
count3 <- count_ss3+count_cc3+count_ECM3

count1 <- con.list[[1]]@net$count[s.cell, s.cell]
count2 <- con.list[[2]]@net$count[s.cell, s.cell]
count3 <- con.list[[3]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2), max(count3))
pdf("Number of interactions_Microglia_End_MOL.pdf", width = 10,height = 5.5)
netVisual_circle(count1, weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12, title.name = paste0("Number of interactions-", names(con.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12, title.name = paste0("Number of interactions-", names(con.list)[2]))
netVisual_circle(count3, weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12, title.name = paste0("Number of interactions-", names(con.list)[3]))
dev.off()

#ligand-receptor
pdf("Compare_signal_bubble_Microglia_MOL_End.pdf", width = 3.5,height =2.5)
levels(cellchat@idents$joint)
netVisual_bubble(cellchat, sources.use = c(7,13), targets.use = c(4), comparison = c(1,2,3), angle.x = 45)
dev.off()

