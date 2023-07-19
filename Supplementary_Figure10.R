
################################################################################
# Supplementary Figure 10
# InN
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

#...TS...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/2.TS')
TS = TS[,TS@meta.data[["celltype"]] %in% c("InN")]
TS$group=str_replace(TS$orig.ident,".*","TS")
save(TS,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/TS.RData")

#...Human...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_1/4.Human')
Human_Fetal = Human[,Human@meta.data[["celltype"]] %in% c("InN","MGE_InN","CGE_InN")]
Human_Fetal$group=str_replace(Human_Fetal$orig.ident,".*","Human_Fetal")
save(Human_Fetal,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Human_Fetal.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_2/3.Human/1.Raw')
Human_Adult1 = Human[,Human@meta.data[["celltype"]] %in% c("MGE_InN","CGE_InN","RC")]
Human_Adult1$group=str_replace(Human_Adult1$orig.ident,".*","Human_Adult1")
save(Human_Adult1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Human_Adult1.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/human/Human_3/0.Raw')
Human_Adult2_DG = Human_DG[,Human_DG@meta.data[["cluster"]] %in% c("InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","InN LHX6 AC008415.1","InN MEIS2 SHISAL2B",
                                                                   "InN NR2F2 ANO2","InN NR2F2 DDR2","InN NR2F2 MIR4300HG","InN NR2F2 PTPRK","InN NR2F2 SLC17A8",
                                                                   "InN PVALB MEPE","InN PVALB PLCL1","InN PVALB PLEKHH2","InN SST ADAMTS12","InN SST EPB41L4A",
                                                                   "InN SST NPY","InN SST OTOF","InN VIP ABI3BP","InN VIP CHRNA2","InN VIP NOX4","InN VIP PENK",
                                                                   "InN VIP SCML4","InN VIP SCTR")]
Human_Adult2_DG$group=str_replace(Human_Adult2_DG$orig.ident,".*","Human_Adult2_DG")
save(Human_Adult2_DG,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Human_Adult2_DG.RData")

Human_CA = Human[,Human@meta.data[["region"]] %in% c("CA1","CA24")]
Human_Adult2_CA = Human_CA[,Human_CA@meta.data[["cluster"]] %in% c("InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","InN LHX6 AC008415.1","InN MEIS2 SHISAL2B",
                                                                   "InN NR2F2 ANO2","InN NR2F2 DDR2","InN NR2F2 MIR4300HG","InN NR2F2 PTPRK","InN NR2F2 SLC17A8",
                                                                   "InN PVALB MEPE","InN PVALB PLEKHH2","InN SST ADAMTS12","InN SST EPB41L4A",
                                                                   "InN SST NPY","InN SST OTOF","InN VIP ABI3BP","InN VIP CHRNA2","InN VIP NOX4","InN VIP PENK",
                                                                   "InN VIP SCML4","InN VIP SCTR")]
Human_Adult2_CA$group=str_replace(Human_Adult2_CA$orig.ident,".*","Human_Adult2_CA")
save(Human_Adult2_CA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Human_Adult2_CA.RData")

#...Macaque...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/monkey/Macaque_3/3.Macaque')
Macaque = Macaque[,Macaque@meta.data[["celltype"]] %in% c("InN")]
Macaque$group=str_replace(Macaque$orig.ident,".*","Macaque")
save(Macaque,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Macaque.RData")


#...Pig...
Pig <- snRNA[,snRNA@meta.data[["age"]] %in% c("Pig")]
Pig1 = Pig[,Pig@meta.data[["seurat_clusters"]] %in% c("2","22","27","32")]
Pig1$group=str_replace(Pig1$orig.ident,".*","Pig1")
save(Pig1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Pig1.RData")

Pig <- sc[,Pig@meta.data[["seurat_clusters"]] %in% c("2")]
Pig = Pig[,Pig@meta.data[["cluster"]] %in% c("InN")]
Pig$group=str_replace(Pig$orig.ident,".*","Pig")
save(Pig,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Pig.RData")

#...Mouse...
setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_3') #弃用
Mouse3 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("5","12")]
Mouse3$group=str_replace(Mouse3$orig.ident,".*","Mouse3")
save(Mouse3,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Mouse3.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_4')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse4 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("4","15","31","35")]
Mouse4$group=str_replace(Mouse4$orig.ident,".*","Mouse4")
Mouse4$orig.ident <- "Mouse_P7"
save(Mouse4,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Mouse4.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/Mouse_2/0.Raw')
Mouse2 = Mouse[,Mouse@meta.data[["characteristics..cell.cluster"]] %in% c("GABA","Immature-GABA")]
Mouse2$group=str_replace(Mouse2$orig.ident,".*","Mouse2")
save(Mouse2,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Mouse2.RData")

setwd('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/mouse/1_Adult')
Mouse <- snRNA[,snRNA@meta.data[["age"]] %in% c("Mouse")]
Mouse1 = Mouse[,Mouse@meta.data[["seurat_clusters"]] %in% c("2","34")]
Mouse1$group=str_replace(Mouse1$orig.ident,".*","Mouse1")
save(Mouse1,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/Mouse1.RData")

##
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Human_Fetal.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Human_Adult1.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Human_Adult2_CA.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Human_Adult2_DG.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Macaque3.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Pig2.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/Mouse3.RData")
load("/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/1.Raw/0.Raw/TS.RData")

################################################################################
# Integrating InN data from different species via using liger method.
################################################################################
snRNA <- merge(TS,y=c(Human_Fetal, Human_Adult1, Human_Adult2_CA,Human_Adult2_DG, Macaque, Mouse1, Mouse4, Pig))
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
save(snRNA,file="/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/InN/cross_species_InN.RData")

################################################################################
# Supplementary Figure 10A
################################################################################
library(Seurat)
table(snRNA$group)  
av<-AverageExpression(snRNA,group.by = "group", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_InN.csv") #保存结果
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))

################################################################################
# Supplementary Figure 10B
################################################################################
markers.to.plot <- c("ERBB4","ADARB2","ZNF385D","GALNTL6","GRIP1","NRXN3")
pdf("TS_Human_MAcaque_special_TS_1.pdf", width = 6,height = 4)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "group", dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# Supplementary Figure 10C
################################################################################
markers.to.plot <- c("ERBB4","ADARB2","ZNF385D","GALNTL6","GRIP1","NRXN3")
pdf("TS_Human_nonHuman_special_TS_3.pdf", width = 6,height = 5)
DotPlot(TS_snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# Supplementary Figure 10D
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

snRNA<- snRNA[,snRNA@meta.data$celltype %in% c("InN")]
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 20, verbose = FALSE)
snRNA <- FindNeighbors(snRNA, dims = 1:10)
snRNA <- FindClusters(snRNA, resolution = 0.1)
snRNA <- RunUMAP(snRNA, dims = 1:20)
snRNA <- RunTSNE(snRNA, dims = 1:20)
head(Idents(snRNA), 5) 
save(snRNA,file="snRNA.RData")
InN_marker<-VlnPlot(snRNA,
                    features = c("CNR1","SST","CALB2","PVALB","RELN","LAMP5","SV2C","CCK","GAD1","GABBR2","SLC6A1","DCX","SNAP25"),group.by = "celltype",stacked=T,pt.size=0,combine = FALSE)+
  theme(strip.text=element_text(size=6),
        axis.title=element_text(size=6,vjust = 1),
        axis.ticks=element_line(size=0),
        axis.text.y=element_text(size=0,colour ='black'),
        axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6),
        axis.title.y=element_text(size=0),
        panel.spacing = unit(-0.1, "lines"))
ggsave(filename = "InN_marker.pdf", plot = cell_marker, device = 'pdf', width = 6, height = 10,units = 'cm')

current.cluster.ids <- c("0","1","2","3","4")
new.cluster.ids <- c("CNR1+SST-","CNR1-SST+","CNR-SST-","RC1","RC2") 
snRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
snRNA <- RenameIdents(snRNA, `0`="CNR1+SST-",`1`="CNR1-SST+",`2`="CNR-SST-",`3`="RC1",`4`="RC2")
Idents(snRNA) <- factor(Idents(snRNA), levels=c("CNR1+SST-","CNR1-SST+","CNR-SST-","RC1","RC2"))

#FindAllMarkers
DefaultAssay(snRNA)="RNA"
snRNA <- NormalizeData(snRNA)
snRNA <- FindVariableFeatures(snRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(snRNA)
snRNA <- ScaleData(snRNA, features = all.genes)
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.01, min.pct = 0.1, only.pos = FALSE, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)

save(snRNA,file="InN.RData")

