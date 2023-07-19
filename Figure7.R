
################################################################################
# Figure 7
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

###...OPC...
################################################################################
# OPC data were extracted for each species.
################################################################################

#...TS...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/2.TS')
table(TS@meta.data[["celltype"]])
TS = TS[,TS@meta.data[["celltype"]] %in% c("OPC")]
TS$group=str_replace(TS$orig.ident,".*","TS")
save(TS,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/TS.RData")

#...Human...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_1/4.Human')
table(Human@meta.data[["celltype"]])
Human_Fetal = Human[,Human@meta.data[["celltype"]] %in% c("OPC")]
Human_Fetal$group=str_replace(Human_Fetal$orig.ident,".*","Human_Fetal")
save(Human_Fetal,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Human_Fetal.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/3.Human/1.Raw')
table(Human@meta.data[["celltype"]])
Human_Adult1 = Human[,Human@meta.data[["celltype"]] %in% c("OPC")]
Human_Adult1$group=str_replace(Human_Adult1$orig.ident,".*","Human_Adult1")
save(Human_Adult1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Human_Adult1.RData")
rm(Human)

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_3/0.Raw')
table(Human_DG@meta.data[["cluster"]])
Human_Adult2_DG = Human_DG[,Human_DG@meta.data[["cluster"]] %in% c("OPC PDGFRA EGR1","OPC PDGFRA GRIA4")]
Human_Adult2_DG$group=str_replace(Human_Adult2_DG$orig.ident,".*","Human_Adult2_DG")
save(Human_Adult2_DG,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Human_Adult2_DG.RData")
rm(Human_DG)

Human_CA = Human[,Human@meta.data[["region"]] %in% c("CA1","CA24")]
table(Human_CA@meta.data[["cluster"]])
Human_Adult2_CA = Human_CA[,Human_CA@meta.data[["cluster"]] %in% c("OPC PDGFRA EGR1","OPC PDGFRA GRIA4")]
Human_Adult2_CA$group=str_replace(Human_Adult2_CA$orig.ident,".*","Human_Adult2_CA")
save(Human_Adult2_CA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Human_Adult2_CA.RData")
rm(Human_CA, Human)

#...Macaque...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/monkey/Macaque_3/3.Macaque')
table(Macaque@meta.data[["celltype"]])
Macaque3 = Macaque[,Macaque@meta.data[["celltype"]] %in% c("OPC")]
Macaque3$group=str_replace(Macaque3$orig.ident,".*","Macaque3")
save(Macaque3,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Macaque3.RData")
#...Macaque...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/monkey/Macaque_2/3.Macaque')
table(Macaque@meta.data[["seurat_clusters"]])
Macaque2 = Macaque[,Macaque@meta.data[["seurat_clusters"]] %in% c("3")]
Macaque2$group=str_replace(Macaque2$orig.ident,".*","Macaque2")
save(Macaque2,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Macaque2.RData")
rm(Macaque)

#...Pig...
Pig <- snRNA[,snRNA@meta.data[["age"]] %in% c("Pig")]  ##华大质量太差
Pig1 = Pig[,Pig@meta.data[["seurat_clusters"]] %in% c("2","22","27","32")]
Pig1$group=str_replace(Pig1$orig.ident,".*","Pig1")
save(Pig1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Pig1.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/pig/Pig_2/3.Pig')
table(Pig@meta.data[["cluster"]])
Pig = Pig[,Pig@meta.data[["cluster"]] %in% c("OPC")]
Pig$group=str_replace(Pig$orig.ident,".*","Pig")
save(Pig,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Pig.RData")

#...Mouse...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_3') #弃用
Mouse = snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse3 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("5","12")]
Mouse3$group=str_replace(Mouse3$orig.ident,".*","Mouse3")
save(Mouse3,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Mouse3.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_4')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse4 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("5","10","12")]
Mouse4$group=str_replace(Mouse4$orig.ident,".*","Mouse4")
Mouse4$orig.ident <- "Mouse_P7"
save(Mouse4,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Mouse4.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_2/0.Raw') ##弃用  太差了
table(Mouse@meta.data[["characteristics..cell.cluster"]])
Mouse2 = Mouse[,Mouse@meta.data[["characteristics..cell.cluster"]] %in% c("OPC")]
Mouse2$group=str_replace(Mouse2$orig.ident,".*","Mouse2")
save(Mouse2,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Mouse2.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/1_Adult')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse1 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("15","24","27")]
Mouse1$group=str_replace(Mouse1$orig.ident,".*","Mouse1")
save(Mouse1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/1.Raw/Mouse1.RData")
rm(Mouse,snRNA)

################################################################################
# Integrating OPC data from different species via using liger method.
################################################################################
snRNA <- merge(TS,y=c(Human_Fetal,Human_Adult1, Human_Adult2_CA, Human_Adult2_DG, Macaque2, Macaque3, Mouse1, Mouse4, Pig))
snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA)
snRNA <- ScaleData(snRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20
snRNA <- RunOptimizeALS(snRNA, k=nFactors, split.by="orig.ident")
snRNA <- RunQuantileNorm(snRNA, split.by="orig.ident")
snRNA$clusters <- factor(snRNA$clusters, levels=1:length(levels(snRNA$clusters)))
snRNA <- FindNeighbors(snRNA, reduction="iNMF", dims=1:nFactors)
snRNA <- FindClusters(snRNA, resolution = 0.5)
snRNA <- RunUMAP(snRNA, dims=1:nFactors,reduction="iNMF")
snRNA <- RunTSNE(snRNA, dims=1:nFactors,reduction="iNMF",check_duplicates = FALSE)
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE, check_duplicates = FALSE)
snRNA_seurat <- ligerToSeurat(snRNA)
DimPlot(seurat_obj, pt.size = 1, label=T, label.size = 4, repel=T)
save(snRNA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/OPC/cross_species_OPC.RData")

################################################################################
# Figure 7A
################################################################################
table(snRNA$group)  
av<-AverageExpression(snRNA,group.by = "group", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_OPC.csv") #保存结果
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))


###...MOL...
################################################################################
# MOL data were extracted for each species.
################################################################################
#...TS...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/2.TS')
TS = TS[,TS@meta.data[["celltype"]] %in% c("MOL")]
TS$group=str_replace(TS$orig.ident,".*","TS")
save(TS,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/TS.RData")

#...Human...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_1/4.Human')
Human_Fetal = Human[,Human@meta.data[["celltype"]] %in% c("InN","MGE_InN","CGE_InN")]
Human_Fetal$group=str_replace(Human_Fetal$orig.ident,".*","Human_Fetal")
save(Human_Fetal,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Human_Fetal.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/3.Human/1.Raw')
Human_Adult1 = Human[,Human@meta.data[["celltype"]] %in% c("NFOL")]
Human_Adult1$group=str_replace(Human_Adult1$orig.ident,".*","Human_Adult1")
save(Human_Adult1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Human_Adult1.RData")
rm(Human)
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_3/0.Raw')
table(Human_DG@meta.data[["cluster"]])
Human_Adult2_DG = Human_DG[,Human_DG@meta.data[["cluster"]] %in% c("Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11")]
Human_Adult2_DG$group=str_replace(Human_Adult2_DG$orig.ident,".*","Human_Adult2_DG")
save(Human_Adult2_DG,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Human_Adult2_DG.RData")
rm(Human_DG)

Human_CA = Human[,Human@meta.data[["region"]] %in% c("CA1","CA24")]
table(Human_CA@meta.data[["cluster"]])
Human_Adult2_CA = Human_CA[,Human_CA@meta.data[["cluster"]] %in% c("Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11")]
Human_Adult2_CA$group=str_replace(Human_Adult2_CA$orig.ident,".*","Human_Adult2_CA")
save(Human_Adult2_CA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Human_Adult2_CA.RData")
rm(Human_CA, Human)

#...Macaque...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/monkey/Macaque_3/3.Macaque')
table(Macaque@meta.data[["celltype"]])
Macaque3 = Macaque[,Macaque@meta.data[["celltype"]] %in% c("NFOL")]
Macaque3$group=str_replace(Macaque3$orig.ident,".*","Macaque3")
save(Macaque3,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Macaque3.RData")
#...Macaque...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/monkey/Macaque_2/3.Macaque')
table(Macaque@meta.data[["seurat_clusters"]])
Macaque2 = Macaque[,Macaque@meta.data[["seurat_clusters"]] %in% c("8","12","29","42")]
Macaque2$group=str_replace(Macaque2$orig.ident,".*","Macaque2")
save(Macaque2,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Macaque2.RData")
rm(Macaque)

#...Pig...
Pig <- snRNA[,snRNA@meta.data[["age"]] %in% c("Pig")]
Pig1 = Pig[,Pig@meta.data[["seurat_clusters"]] %in% c("2","22","27","32")]
Pig1$group=str_replace(Pig1$orig.ident,".*","Pig1")
save(Pig1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Pig1.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/pig/Pig_2/3.Pig')
table(Pig@meta.data[["clusters"]])
Pig = Pig[,Pig@meta.data[["cluster"]] %in% c("Oligo")]
Pig$group=str_replace(Pig$orig.ident,".*","Pig")
save(Pig,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Pig.RData")

#...Mouse...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_3') #弃用
Mouse = snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse3 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("5","12")]
Mouse3$group=str_replace(Mouse3$orig.ident,".*","Mouse3")
save(Mouse3,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Mouse3.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_4')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse4 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("6","14","16","33","34")]
Mouse4$group=str_replace(Mouse4$orig.ident,".*","Mouse4")
Mouse4$orig.ident <- "Mouse_P7"
save(Mouse4,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Mouse4.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_2/0.Raw')
table(Mouse@meta.data[["characteristics..cell.cluster"]])
Mouse2 = Mouse[,Mouse@meta.data[["characteristics..cell.cluster"]] %in% c("MOL")]
Mouse2$group=str_replace(Mouse2$orig.ident,".*","Mouse2")
save(Mouse2,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Mouse2.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/1_Adult')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse1 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("4","9")]
Mouse1$group=str_replace(Mouse1$orig.ident,".*","Mouse1")
save(Mouse1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/1.Raw/Mouse1.RData")

################################################################################
# Integrating MOL data from different species via using liger method.
################################################################################
snRNA <- merge(TS,y=c(Human_Adult1, Human_Adult2_CA, Human_Adult2_DG, Macaque2, Macaque3, Mouse1, Mouse4, Pig))
snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA)
snRNA <- ScaleData(snRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20 
snRNA <- RunOptimizeALS(snRNA, k=nFactors, split.by="orig.ident")
snRNA <- RunQuantileNorm(snRNA, split.by="orig.ident")
snRNA$clusters <- factor(snRNA$clusters, levels=1:length(levels(snRNA$clusters)))
snRNA <- FindNeighbors(snRNA, reduction="iNMF", dims=1:nFactors)
snRNA <- FindClusters(snRNA, resolution = 0.4)
snRNA <- RunUMAP(snRNA, dims=1:nFactors,reduction="iNMF")
snRNA <- RunTSNE(snRNA, dims=1:nFactors,reduction="iNMF",check_duplicates = FALSE)
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE, check_duplicates = FALSE)
snRNA_seurat <- ligerToSeurat(snRNA)
DimPlot(seurat_obj, pt.size = 1, label=T, label.size = 4, repel=T)
save(snRNA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/MOL/cross_species_MOL.RData")

################################################################################
# Figure 7B
################################################################################
table(snRNA$group)  
av<-AverageExpression(snRNA,group.by = "group", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_MOL.csv") #保存结果
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))

################################################################################
# Figure 7C
################################################################################
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)

load("~/Downloads/snRNA/snRNA.RData")
snRNA <- snRNA[,snRNA@meta.data$celltype %in% c("OPC")]
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
count_table <- table(snRNA@meta.data[["seurat_clusters"]], snRNA@meta.data[["Age"]])
count_table
current.cluster.ids <- c("0","1")
new.cluster.ids <- c("OPC1","OPC2") 
snRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Idents(snRNA) <- snRNA$celltype
Idents(snRNA) <- factor(Idents(snRNA), levels=c("OPC1","OPC2"))
save(snRNA,file="TS_OPC.RData")
count_table <- table(snRNA@meta.data[["celltype"]], snRNA@meta.data[["Age"]])
count_table
p1 <- DimPlot(snRNA, reduction = "umap", group.by = "Age", pt.size=1)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p2 <- DimPlot(snRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap1 <- plot_grid(p1, p2,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 24, height = 9, units = 'cm')
#
DefaultAssay(snRNA) <- "RNA"
OPC_marker <- VlnPlot(snRNA,features = c("OLIG2","OLIG1","PDGFRA","VCAN"), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6))
ggsave(filename = "OPC_marker.pdf", plot = OPC_marker, device = 'pdf', width = 5, height = 4, units = 'cm')
rm('OPC_marker')

################################################################################
# Figure 7E
################################################################################
DefaultAssay(snRNA)="RNA"
snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(snRNA)
snRNA <- ScaleData(snRNA, features = all.genes)
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "MAST")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
DefaultAssay(snRNA) <- "RNA"
markers.to.plot <- c("CNTN6","IGF1R","NINJ2","TTR","TMEM98","DDB2","SNTG1","OLIG1","OLIG2","DDIT4","RALYL")
pdf("OPC_DEGs.pdf", width = 6,height = 1.7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# Figure 7F
################################################################################
snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("MOL")]
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE)
ElbowPlot(snRNA)
snRNA <- FindNeighbors(snRNA, dims = 1:5)
snRNA <- FindClusters(snRNA, resolution = 0.2)
snRNA <- RunUMAP(snRNA, dims = 1:10)
snRNA <- RunTSNE(snRNA, dims = 1:10) 
head(Idents(snRNA), 5) 
current.cluster.ids <- c("0","1","2")
new.cluster.ids <- c("MOL1","MOL2","MOL3") 
snRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
snRNA <- RenameIdents(snRNA, `0`="MOL1",`1`="MOL2",`2`="MOL3")
Idents(snRNA) <- factor(Idents(snRNA), levels=c("MOL1","MOL2","MOL3"))
save(snRNA,file="TS_MOL.RData")
count_table <- table(snRNA@meta.data[["celltype"]], snRNA@meta.data[["Age"]])
count_table
p1 <- DimPlot(snRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = F,repel = F)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p2 <- DimPlot(snRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
umap1 <- plot_grid(p1, p2,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 23.5, height = 8, units = 'cm')
DefaultAssay(snRNA) <- "RNA"
MOL_marker <- VlnPlot(snRNA,features = c("MBP","MOG","PLP1","ST18"), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6))
ggsave(filename = "MOL_marker.pdf", plot = MOL_marker, device = 'pdf', width = 5, height = 4, units = 'cm')
rm('MOL_marker')

################################################################################
# Figure 7H
################################################################################
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "MAST")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
DefaultAssay(snRNA) <- "RNA"
markers.to.plot <- c("CNTN6","IGF1R","NINJ2","TTR","TMEM98","DDB2","SNTG1","OLIG1","OLIG2","DDIT4","RALYL")
pdf("MOL_DEGs.pdf", width = 6,height = 1.7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# Figure 7J-N
################################################################################
# monocle2
rm(list=ls())
library(Seurat)；library(monocle)
load("~/Downloads/snRNA/snRNA.RData")
snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("OPC","NFOL","MOL")]
snRNA <- subset(snRNA, downsample = 2000)
data <- as(as.matrix(snRNA@assays[["RNA"]]@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = snRNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
cds  <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,expressionFamily = negbinomial.size())
HSMM <- cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 1 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 1))
print(head(pData(HSMM)))
#Clustering cells without marker genes 
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.001)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.05)) 
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)
HSMM@auxOrderingData[[HSMM@dim_reduce_type]]$branch_points
count_table <- table(HSMM$age, HSMM$State)
count_table

################################################################################
# Figure 7J
################################################################################
plot_cell_trajectory(HSMM,color_by = "celltype",cell_size = 0.1)
ggsave("celltype.pdf",device = "pdf",width = 10,height = 12,units = c("cm"))
plot_cell_trajectory(HSMM, color_by = "State",cell_size = 0.2)+facet_wrap(~State, nrow = 1)
ggsave("state.pdf",device = "pdf",width = 21,height = 9,units = c("cm"))

################################################################################
# Figure 7K
################################################################################
HSMM_genes <- row.names(subset(fData(HSMM),gene_short_name %in% c("PDGFRA","ST18")))
plot_genes_branched_pseudotime(HSMM[HSMM_genes,],cell_size = 0.01,color_by = "Pseudotime",ncol = 1)
ggsave("OPC_MOL_pseudotime.pdf",device = "pdf",width = 10,height = 9,units = c("cm"))
plot_genes_branched_pseudotime(HSMM[HSMM_genes,],cell_size=0.01, color_by = "State",ncol = 1)
ggsave("OPC_MOL_State.pdf",device = "pdf",width = 10,height = 9,units = c("cm"))
plot_genes_branched_pseudotime(HSMM[HSMM_genes,],cell_size=0.01, color_by = "celltype",ncol = 1)
ggsave("OPC_MOL_celltype.pdf",device = "pdf",width = 10,height = 9,units = c("cm"))

################################################################################
# Figure 7L
################################################################################
SMM = orderCells(HSMM,root_state = 2) 
BEAM_res = BEAM(HSMM,branch_point = 1,cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
saveRDS(BEAM_res, file = "BEAM_res1.rds")
library('RColorBrewer')
tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 3, #这些基因被分成几个group
                                 cores = 1,
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T
)
pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()

#GO
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
write.csv(gene_group,"GP_gene_group.csv") #保存结果
write.csv(allcluster_go,"allcluster_go.csv") #保存结果
write.csv(go_res,"go_res.csv") #保存结果
write.csv(diff_test_res,"diff_test_res.csv") #保存结果
saveRDS(HSMM, file = "HSMM.rds")

################################################################################
# Figure 7M
################################################################################
plot_complex_cell_trajectory(HSMM, x=1, y=2, color_by = "State",cell_size = 0.01)+theme(legend.title = element_blank())+facet_wrap(~age,nrow=1) #树形图
ggsave("tree.pdf",device = "pdf",width = 12,height = 6.5,units = c("cm"))

################################################################################
# Figure 7N
################################################################################
HSMM_genes <- row.names(subset(fData(HSMM),gene_short_name %in% c("CNTN6","IGF1R","NINJ2","DDB2","TMEM98","IL33","SNTG1","RALYL","SOX2","HMGN3","HAPLN2","TNFAIP6")))
plot_genes_branched_pseudotime(HSMM[HSMM_genes,],cell_size=0.01, color_by = "State",ncol = 1)
ggsave("OPC_MOL_aging_State.pdf",device = "pdf",width = 10,height = 45,units = c("cm"))

