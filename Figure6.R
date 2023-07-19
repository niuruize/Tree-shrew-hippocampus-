
###...Microglia...
################################################################################
# Figure 6A
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
# Microglia data were extracted for each species.
################################################################################
#...TS...
TS_Microglia = TS[,TS@meta.data[["celltype"]] %in% c("Microglia")]
TS_Microglia$group=str_replace(TS_Microglia$orig.ident,".*","TS_Microglia")
save(TS_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/TS_Microglia.RData")

#...Human...
Human_Fetal_Microglia = Human[,Human@meta.data[["celltype"]] %in% c("Microglia")]
Human_Fetal_Microglia$group=str_replace(Human_Fetal_Microglia$orig.ident,".*","Human_Fetal_Microglia")
save(Human_Fetal_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Human_Fetal_Microglia.RData")
Human_Adult1_Microglia = Human[,Human@meta.data[["celltype"]] %in% c("Microglia")]
Human_Adult1_Microglia$group=str_replace(Human_Adult1_Microglia$orig.ident,".*","Human_Adult1_Microglia")
save(Human_Adult1_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Human_Adult1_Microglia.RData")

Human_Adult2_DG_Microglia = Human_DG[,Human_DG@meta.data[["cluster"]] %in% c("Micro C1QB CD83","Micro C1QB P2RY12")]
Human_Adult2_DG_Microglia$group=str_replace(Human_Adult2_DG_Microglia$orig.ident,".*","Human_Adult2_DG_Microglia")
save(Human_Adult2_DG_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Human_Adult2_DG_Microglia.RData")

Human_Adult2_CA_Microglia = Human_CA[,Human_CA@meta.data[["cluster"]] %in% c("Micro C1QB CD83","Micro C1QB P2RY12")]
Human_Adult2_CA_Microglia$group=str_replace(Human_Adult2_CA_Microglia$orig.ident,".*","Human_Adult2_CA_Microglia")
save(Human_Adult2_CA_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Human_Adult2_CA_Microglia.RData")

#...Macaque...
Macaque_Microglia = Macaque[,Macaque@meta.data[["celltype"]] %in% c("Microglia")]
Macaque_Microglia$group=str_replace(Macaque_Microglia$orig.ident,".*","Macaque_Microglia")
save(Macaque_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Macaque_Microglia.RData")

#...Pig...
Pig <- snRNA[,snRNA@meta.data[["age"]] %in% c("Pig")]
Pig_Microglia = Pig[,Pig@meta.data[["seurat_clusters"]] %in% c("2","22","27","32")]
Pig_Microglia$group=str_replace(Pig_Microglia$orig.ident,".*","Pig_Microglia")
save(Pig_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Pig_Microglia.RData")

#...Mouse...
Mouse1_Microglia = Mouse1[,Mouse1@meta.data[["seurat_clusters"]] %in% c("5","12")]
Mouse1_Microglia$group=str_replace(Mouse1_Microglia$orig.ident,".*","Mouse1_Microglia")
save(Mouse1_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Mouse1_Microglia.RData")

Mouse2_Microglia = Mouse2[,Mouse2@meta.data[["cluster_name"]] %in% c("Microglia")]
Mouse2_Microglia$group=str_replace(Mouse2_Microglia$orig.ident,".*","Mouse2_Microglia")
save(Mouse2_Microglia,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/Mouse2_Microglia.RData")

################################################################################
# Integrating microglial data from different species via using liger method.
################################################################################
snRNA <- merge(TS_Microglia,y=c(Human_Fetal_Microglia, Human_Adult1_Microglia, Human_Adult2_CA_Microglia,Human_Adult2_DG_Microglia,
                                Macaque_Microglia, Mouse1_Microglia, Mouse2_Microglia, Pig_Microglia))

snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA)
snRNA <- ScaleData(snRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20 
snRNA <- RunOptimizeALS(snRNA, k=nFactors, split.by="orig.ident")
snRNA <- RunQuantileNorm(snRNA, split.by="orig.ident")
snRNA$clusters <- factor(snRNA$clusters, 
                         levels=1:length(levels(snRNA$clusters)))
snRNA <- FindNeighbors(snRNA, reduction="iNMF", dims=1:nFactors)
snRNA <- FindClusters(snRNA, resolution = 0.4)
snRNA <- RunUMAP(snRNA, dims=1:nFactors,reduction="iNMF")
snRNA <- RunTSNE(snRNA, dims=1:nFactors,reduction="iNMF")
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE)
snRNA_seurat <- ligerToSeurat(snRNA)
DimPlot(seurat_obj, pt.size = 1, label=T, label.size = 4, repel=T)
save(snRNA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/cross_species_Micro.RData")

################################################################################
# Figure 6A
################################################################################
DimPlot(snRNA, reduction = "umap",group.by = "group", 
        cols =c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282'), pt.size=0.1)+
  theme(plot.title=element_text(size=0),plot.tag =element_text(size=0), strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_group.pdf", device = 'pdf', width = 20, height = 14,  units = 'cm')

################################################################################
# Figure 6B
################################################################################
table(snRNA$group)  
av<-AverageExpression(snRNA,group.by = "group", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_microglia.csv")
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))

################################################################################
# Figure 6C
################################################################################
Microglia_specific_marker<-VlnPlot(snRNA,group.by = "group",features = c("CD33","GRN","TYROBP","TREM2","APOE","CD68","CD81","C1QC","C1QB","CTSS","CX3CR1","LAPTM5","P2RY12","PTPRC","CSF1R","MEF2C"), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=90,size=6))
ggsave(filename = "Microglia_specific_marker.pdf", plot = Microglia_specific_marker, device = 'pdf', width = 12, height = 20, units = 'cm')
rm('all_marker')

################################################################################
# Figure 6D
################################################################################
# SCENIC
# SCENIC workflow:
load('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/cross_species_Micro.RData')
library(Seurat)
AB1 <- subset(snRNA, downsample = 1000)
table(AB1@meta.data[["group"]])
exprMat<-GetAssayData(object = AB1)#count
dim(exprMat)
exprMat[1:4,1:10]
save(exprMat,file = "int/exprMat.Rds")
cellInfo <-data.frame(AB1@meta.data)
head(cellInfo)
groupColumn <- "group"
colnames(cellInfo)[which(colnames(cellInfo)==groupColumn)] <- "group"
cbind(table(cellInfo$group))
head(cellInfo)
saveRDS(cellInfo, file="int/cellInfo.Rds")
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(group=c("Human_Adult1_Microglia"="#E5D2DD",
                        "Human_Adult2_CA_Microglia"="#53A85F",
                        "Human_Adult2_DG_Microglia"="#F1BB72",
                        "Human_Fetal_Microglia"="#F3B1A0",
                        "Macaque_Microglia"="#D6E7A3",
                        "Mouse1_Microglia"="#57C3F3",
                        "Mouse2_Microglia"="#476D87",
                        "Pig_Microglia"="#E95C59",
                        "TS_Microglia"="#AB3282"))
colVars$group <- colVars$group[intersect(names(colVars$group), cellInfo$group)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$group, legend=names(colVars$group))
#
library(SCENIC)
org="hgnc" # or mgi-mouse, or dmel-fly
dbDir="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/5.SCENIC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="Microglia" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
#Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
#Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
options("mc.cores"=10)
#####
library(RcisTarget)
dbFilePath <- getDatabases(scenicOptions)[[1]]
motifRankings <- importRankings(dbFilePath)
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells<-rownames(exprMat)
length(genesLeft_minCells)
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)
genesKept <- genesLeft_minCells_inDatabases
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
class(exprMat_filtered)
runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
library(GENIE3)
runGenie3(as.matrix(exprMat_filtered), scenicOptions)
save(exprMat_filtered,scenicOptions,file = "input_GENIE3_data.Rdata")
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, 
                           minCountsPerGene=3*.01*ncol(exprMat), 
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
save(exprMat,file = "int/exprMat.Rds")
library(GENIE3)
logMat <- log2(exprMat_filtered+1)
runGenie3(logMat, scenicOptions)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget"))
runSCENIC_3_scoreCells(scenicOptions, logMat)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)
#
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType1 <- sapply(split(rownames(cellInfo), cellInfo$group),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType1_Scaled <- t(scale(t(regulonActivity_byCellType1), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType1_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
topRegulators <- reshape2::melt(regulonActivity_byCellType1_Scaled)
colnames(topRegulators) <- c("Regulon", "group", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_cell.csv") #保存结果

minPerc <- .1
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType1_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$group),
                                                function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType1_Binarized[which(rowSums(regulonActivity_byCellType1_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)

# After obtaining the topRegulators data, use cytoscape to build the network graph

################################################################################
# Figure 6L
################################################################################
load("~/Downloads/snRNA/snRNA.RData")
#HuMi_score
HuMi_gene <- read_xlsx("HuMi.xlsx")
gene <- as.list(HuMi_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "HuMi")
colnames(AB@meta.data)[11]<-"HuMi_Score" 
##RidgePlot
P1 <- RidgePlot(AB, features = 'HuMi_Score', ncol = 1) 
ggsave(filename = "HuMi_RidgePlot.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1')  

################################################################################
# Figure 6M
################################################################################
# Scores of HuMi gene sets for all cell type in different groups
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggboxplot(b, x = "orig.ident", y = "HuMi_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "orig.ident", palette = "npg")+stat_compare_means(comparisons = my_comparisons)
ggsave("HuMi_8.pdf",width = 15,height = 15,units = "cm")

# Scores of HuMi gene sets for microglia in different groups
tmp= AB[,AB@meta.data$celltype %in% c("Microglia")]
b<-FetchData(tmp, vars = c("orig.ident","HuMi_Score","celltype"))
b[["orig.ident"]]<-factor(b[["orig.ident"]], levels=c("Infancy","Adult","Old"))
ggviolin(b, x = "orig.ident", y = "HuMi_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("HuMi_Microglia_3.pdf",width = 10,height = 10.5,units = "cm")

################################################################################
# Figure 6N,O
################################################################################
load("~/Downloads/snRNA/snRNA.RData")
snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("Micro")]
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE)
snRNA <- FindNeighbors(snRNA, dims = 1:10)
snRNA <- FindClusters(snRNA, resolution = 0.2)
snRNA <- RunUMAP(snRNA, dims = 1:10)
snRNA <- RunTSNE(snRNA, dims = 1:10) 
head(Idents(snRNA), 5) 
current.cluster.ids <- c("0","1","2")
new.cluster.ids <- c("Micro1","Micro2","Micro3") 
snRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Idents(snRNA) <- snRNA@meta.data$celltype
Idents(snRNA) <- factor(Idents(snRNA), levels=c("Micro1","Micro2","Micro3"))
save(snRNA,file="TS_Micro.RData")
count_table <- table(snRNA@meta.data[["seurat_clusters"]], snRNA@meta.data[["Age"]])
count_table
# plot the count_table using Prism 9
p <- DimPlot(snRNA, reduction = "umap",group.by = "celltype",cols =c('#476D87', '#53A85F', '#F1BB72'), pt.size=0.5)+
  theme(plot.title=element_text(size=0),plot.tag =element_text(size=0), strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "sub_micro_umap.pdf", plot = p, device = 'pdf', width = 11.5, height = 9, units = 'cm')

# FindAllMarkers
DefaultAssay(snRNA)="RNA"
snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(snRNA)
snRNA <- ScaleData(snRNA, features = all.genes)
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)

################################################################################
# Figure 6P,Q
################################################################################
load("~/Downloads/snRNA/TS_Micro.RData")
snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("Micro1")]
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE)
ElbowPlot(snRNA)
snRNA <- FindNeighbors(snRNA, dims = 1:10)
snRNA <- FindClusters(snRNA, resolution = 0.1)
snRNA <- RunUMAP(snRNA, dims = 1:10)
snRNA <- RunTSNE(snRNA, dims = 1:10)
head(Idents(snRNA), 5) 
current.cluster.ids <- c("0","1")
new.cluster.ids <- c("Micro1_1","Micro1_2") 
snRNA$subcelltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Idents(snRNA) <- snRNA$subcelltype
Idents(snRNA) <- factor(Idents(snRNA), levels=c("Micro1_1","Micro1_2"))
count_table <- table(snRNA@meta.data[["subcelltype"]], snRNA@meta.data[["age"]])
count_table
p <- DimPlot(snRNA, reduction = "umap",group.by = "subcelltype", 
             cols =c('#476D87', '#53A85F'), pt.size=0.5)+
  theme(plot.title=element_text(size=0),plot.tag =element_text(size=0), strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "Micro1_subtype_umap.pdf", plot = p11, device = 'pdf', width = 11, height = 9, units = 'cm')
p <- Nebulosa::plot_density(snRNA, features = c("C1QB","C1QA","CX3CR1"),joint = T,reduction = "umap",size = 1)[[4]]
ggsave(filename = "Micro1_2.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(snRNA, features = c("BCL2L1","FYN","GAB2","MAPK1"),joint = T,reduction = "umap",size = 1)[[5]]
ggsave(filename = "Micro1_1.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
# FindAllMarkers
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)

