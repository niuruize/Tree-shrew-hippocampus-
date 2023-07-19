#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

#loading data
load("~/AB.RData")

################################################################################
# Supplementary Figure 1A 
################################################################################

p1 <- DimPlot(AB, reduction = "umap", group.by = "seurat_clusters", pt.size=0.1, label =F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p2 <- DimPlot(AB, reduction = "umap", group.by = "seurat_clusters", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
umap2 <- plot_grid(p1, p1, align = "v",ncol = 2)
ggsave(filename = "seurat_clusters_umap.pdf", plot = umap2, device = 'pdf', width = 31, height = 8, units = 'cm')
p1 <- DimPlot(AB, reduction = "umap", group.by = "celltype", split.by = "Age",pt.size=0.1, label =F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p2 <- DimPlot(AB, reduction = "umap", group.by = "celltype", split.by = "Age",pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
umap2 <- plot_grid(p1, p1, align = "v",ncol = 2)
ggsave(filename = "celltype_age_umap.pdf", plot = umap2, device = 'pdf', width = 31, height = 8, units = 'cm')

################################################################################
###---Supplementary Figure 1B;umap for marker gene---###
################################################################################

celltype_umap<- FeaturePlot(AB, features = c("AQP4","GFAP"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Astro_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("CSF1R","PTPRC"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Micro_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("C1QA","P2RY12"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Micro_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("C1QB","MEF2C"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Micro_umap_3.pdf", plot = celltype_umap, device = 'pdf', width =9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("SLC1A3","SOX2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NPC_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("GFAP",  "SOX6"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NPC_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("PAX6",  "HOPX"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NPC_umap_3.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("VCAN","PDGFRA"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "OPC_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("BCAS1","FYN"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NFOL_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("SNAP25","DCX"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "ImN_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("CAMK2A","PROX1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "DG_ExN_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("CAMK2A","SATB2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "DG_ExN_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("SNAP25","SLC6A1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "InN_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("CCK","LAMP5"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "DG_RC_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("SV2C","CNR1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "DG_RC_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("SST","CALB2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "DG_RC_umap_3.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("RELN","PVALB"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "CR_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("EBF1","RGS5"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "End_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("MOG","ST18"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "MOL_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("DNAH9","CFAP54"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Ependymal_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(AB, features = c("DCN","COL1A2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Pericyte_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)

################################################################################
###---Supplementary Figure 1C; dot plot for marker gene---###
################################################################################

markers.to.plot <- c("SLC1A3","SOX2","GFAP","PAX6","HOPX","AQP4","SOX6","PTPRC","C1QA","C1QB","CSF1R","VCAN","OLIG1","PDGFRA",
                     "BCAS1","FYN","MOG","PLP1","ST18","SNAP25","PROX1","DCX","STMN2","CAMK2A","SLC6A1","CCK","LAMP5","SV2C","CNR1",
                     "SST","CALB2","PVALB","RELN","EBF1","RGS5","DNAH9","CFAP54","DCN","COL1A2")
pdf("celltype_marker.pdf", width = 15.5,height = 7)
DotPlot(AB, features = markers.to.plot, assay='RNA' )+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+coord_flip()
dev.off()

################################################################################
###---Supplementary Figure 1D---###
################################################################################
table(AB$celltype)
# plot the bar graph using Prism 9

################################################################################
###---Supplementary Figure 1E---###
################################################################################
markers <- FindAllMarkers(AB, logfc.threshold = 0, min.pct = 0.1, only.pos = FALSE, test.use = "wilcox")
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
#GSEA
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(enrichplot)
#load data
df = read.table("markers.txt",header = T)
df = df[2:5679,]#NSC
df = df[5680:9060,]#IPC1
df = df[9061:11474,]#IPC2
df = df[11475:17004,]#Microglia
df = df[17005:22386,]#OPC
df = df[22387:25990,]#NFOL
df = df[25991:31791,]#MOL
df = df[31792:38944,]#ImN
df = df[38945:44884,]#DG_ExN
df = df[44885:51746,]#nonDG_ExN
df = df[51747:57973,]#InN
df = df[57974:60940,]#RC
df = df[60941:66542,]#End
df = df[66543:70658,]#Ependymal
df = df[70659:74786,]#Pericyte
df = df[74787:80173,]#Unk
head(df)
dim(df)
#gene SYMBOL to ENTREZID
df_id<-bitr(df$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_all <- merge(df,df_id,by="SYMBOL",all=F)
head(df_all)
dim(df_all)
#GSEA
df_all_sort <- df_all[order(df_all$avg_log2FC, decreasing = T),]
gene_fc = df_all_sort$avg_log2FC
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID
head(gene_fc)
GO <- gseGO(
  gene_fc, 
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
head(GO)
#output
sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]
head(sortGO)
dim(sortGO)
write.csv(sortGO,"NSC.csv")
# plot the heatmap using TBtools

################################################################################
###---Supplementary Figure 1F---###
################################################################################

library(Seurat);library(RcisTarget);library(SCENIC);library(GENIE3)
#load seurat object 
load('/AB.RData')
AB1 <- subset(AB, downsample = 1000)
table(AB1@meta.data[["celltype"]])
exprMat<-GetAssayData(object = AB1)#count
dim(exprMat)
exprMat[1:4,1:10]
save(exprMat,file = "int/exprMat.Rds")
cellInfo <-data.frame(AB1@meta.data)
head(cellInfo)
celltypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==celltypeColumn)] <- "celltype"
cbind(table(cellInfo$celltype))
head(cellInfo)
saveRDS(cellInfo, file="int/cellInfo.Rds")
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(celltype=c("NSC/Astrocyte"="#E5D2DD",
                           "IPC1"="#53A85F",
                           "IPC2"="#F1BB72",
                           "Microglia"="#F3B1A0",
                           "OPC"="#D6E7A3",
                           "NFOL"="#57C3F3",
                           "MOL"="#476D87",
                           "ImN"="#E95C59",
                           "DG_ExN"="#E59CC4",
                           "nonDG_ExN"="#AB3282",
                           "InN"="#23452F",
                           "Ependymal"="#BD956A",
                           "Pericyte"="#8C549C",
                           "CR"="#585658",
                           "End"="#9FA3A8",
                           "Unk"="#E0D4CA"))
colVars$cellcype <- colVars$celltype[intersect(names(colVars$celltype), cellInfo$celltype)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$celltype, legend=names(colVars$celltype))
#
org="hgnc" # or mgi-mouse, or dmel-fly
dbDir="~/SCENIC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="TS" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
#Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
#Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#data filtering
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
save(exprMat,file = "int/exprMat.Rds")
logMat <- log2(exprMat_filtered+1)
logMat <- na.omit(logMat)
save(logMat,file = "int/logMat.Rds")
runGenie3(logMat, scenicOptions)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget")) #** Only for toy run!!
runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, logMat)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType1 <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType1_Scaled <- t(scale(t(regulonActivity_byCellType1), center = T, scale=T))
pdf("celltype_TF_2.pdf", width = 6,height = 9)
pheatmap::pheatmap(regulonActivity_byCellType1_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
topRegulators <- reshape2::melt(regulonActivity_byCellType1_Scaled)
colnames(topRegulators) <- c("Regulon", "celltype", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_cell.csv")
# plot the network using cytoscape















