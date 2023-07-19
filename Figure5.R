################################################################################
#---Figure 5
################################################################################
#loading R packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)

load("~/Downloads/snRNA/cross_species_snRNA.RData")
snRNA$celltype_species <- paste(snRNA$celltype, snRNA$species, sep = "_")
table(snRNA$celltype_species)
snRNA<- snRNA[,scRNA@meta.data$celltype_species %in% c("Human_NSC/Astrocyte","Macaque_NSC/Astrocyte","Mouse_Astrocyte","Mouse_nIPC","Mouse_nIPC_p","Mouse_RGL","TS_IPC1","TS_IPC2","TS_NSC/Astrocyte")]
table(snRNA$orig.ident)
snRNAlist <- SplitObject(snRNA, split.by = "orig.ident")
snRNAlist <- lapply(X = snRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
snRNA.features <- SelectIntegrationFeatures(object.list = snRNAlist, nfeatures = 3000)
snRNAlist <- lapply(X = snRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = snRNA.features, verbose = FALSE)
  x <- RunPCA(x, features = snRNA.features, verbose = FALSE)
})
snRNA.anchors <- FindIntegrationAnchors(object.list = snRNAlist, normalization.method = "LogNormalize",reduction = "rpca", anchor.features = snRNA.features, dims = 1:50)  ##耗时久
snRNA <- IntegrateData(snRNA.anchors, normalization.method = "LogNormalize", dims = 1:50)  ##耗时久
snRNA
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 30, verbose = FALSE)
snRNA <- JackStraw(snRNA, num.replicate = 100)  ##耗时久
snRNA <- ScoreJackStraw(snRNA, dims = 1:20) 
JackStrawPlot(snRNA, dims = 1:20)
ElbowPlot(snRNA)
snRNA <- FindNeighbors(snRNA, dims = 1:10)
snRNA <- FindClusters(snRNA, resolution = 0.4)
snRNA <- RunUMAP(snRNA, dims = 1:10)
snRNA <- RunTSNE(snRNA, dims = 1:10)
head(Idents(snRNA), 5)
markers.to.plot <- c("SLC1A3","GFAP","AQP4","HES5","SOX2","EMX2","FOXG1","EGFR","LPAR1","ETNPPL",
                     "MFGE8","SOX5","SOX6","CDK1","ASCL1","TFAP2C","EOMES","IGFBPL1","CALB2","PLK5",
                     "NES","NEUROG2","MCM2","DCX","NEUROD1","CALB1","TOP2A","ELAVL2","VIM","NOTCH2",
                     "SOX4","SOX11", "PADI2","FRZB","WNT8B","MAP2","ADAMTS19","CDK6","WIF1")
pdf("celltype_marker.pdf", width = 15.5,height = 7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
#
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")
new.cluster.ids <- c("Astrocyte","Astrocyte","Astrocyte","RGL1","RGL1","nIPC_p","Astrocyte","NB","RGL2","TS_RGL","Immu_RGL") 
snRNA@meta.data$subcelltype <- plyr::mapvalues(x = as.integer(as.character(snRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
snRNA@meta.data[["subcelltype"]] <- factor(snRNA$subcelltype, levels=c("Astrocyte","RGL1","RGL2","TS_RGL","nIPC_p","NB","Immu_RGL"))
snRNA@meta.data[["species"]]<-factor(snRNA@meta.data[["species"]], levels=c("Human","Macaque","TS","Mouse"))
save(snRNA,file="cross_species_NSC.RData")

################################################################################
#---Figure 5A
################################################################################
p <- DimPlot(snRNA, reduction = "umap", group.by = "group", pt.size=0.01)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "species_umap.pdf", plot = p, device = 'pdf', width = 24, height = 9, units = 'cm')

################################################################################
#---Figure 5B
################################################################################
p <- DimPlot(snRNA, reduction = "umap", group.by = "subcelltype", pt.size=0.01)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "subcelltype_umap.pdf", plot = p, device = 'pdf', width = 24, height = 9, units = 'cm')

###############################################################################
#---Figure 5C
################################################################################
snRNA$subcelltype_species <- paste(snRNA$subcelltype, snRNA$species, sep = "_")
p <- DimPlot(snRNA, reduction = "umap", group.by = "subcelltype_species", pt.size=0.01)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "subcelltype_species_umap.pdf", plot = p, device = 'pdf', width = 24, height = 9, units = 'cm')

###############################################################################
#---Figure 5D
################################################################################
DefaultAssay(snRNA) <- "RNA"
celltype_umap<- FeaturePlot(snRNA, features = c("AQP4","GFAP"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Astro_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("MFGE8","ASCL1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("VIM","PDGFRB"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("NES","PROX1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_3.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("SLC1A3","SOX2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_4.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("CDK1","TOP2A"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_5.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("NEUROD1","CALB2"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_6.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("SOX4","SOX11"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_7.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("SOX5","SOX6"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_8.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("MAP2","ADAMTS19"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_9.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("CDK6","WIF1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_10.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(snRNA, features = c("DCX","IGFBPL1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NB_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)

###############################################################################
#---Figure 5E
################################################################################
DefaultAssay(snRNA) <- "RNA"
p <- Nebulosa::plot_density(snRNA, features = c("ADAMTS19","MAP2","CDK6","WIF1"),joint = T,reduction = "umap",size = 0.01)[[5]]
ggsave(filename = "TS_RGL.pdf", plot = p, device = 'pdf', width = 10, height = 7.5, units = 'cm')

###############################################################################
#---Figure 5F
################################################################################
Idents(snRNA)="subcelltype"
markers <- FindAllMarkers(snRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE, test.use = "wilcox")
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
top30 = markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
TS_RGL_top30 = subset(top30, top30$cluster==TS_RGL)
pdf("TS_RGL_top30_DEGs.pdf", width = 15.5,height = 7)
DotPlot(snRNA, features = TS_RGL_top30$gene, cols = c("blue", "red","green"), group.by = "subcelltype", dot.scale = 8) +  RotatedAxis()
dev.off()

