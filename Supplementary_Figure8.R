################################################################################
# Supplementary Figure 8 C
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

p1<-VlnPlot(AB, features = c("CREBBP"),idents = 'NSC/Astrocyte',pt.size = 0.1,split.by = "age",ncol = 1)+theme(axis.title.x=element_text(size=0,angle=0))
ggsave(filename = "NSC_Astrocyte_SLC1A3.pdf", plot = p1, device = 'pdf', width = 5, height = 5)

p1<-VlnPlot(AB, features = c("HDAC9"),idents = 'NSC/Astrocyte', pt.size = 0.1,split.by = "age",ncol = 1)+theme(axis.title.x=element_text(size=0,angle=0))
ggsave(filename = "DG_ExN_SNAP25.pdf", plot = p1, device = 'pdf', width = 5, height = 5)

################################################################################
# Supplementary Figure 8 G
################################################################################
pdf("Astrocyte_DEGs.pdf", width = 5,height = 2)
DotPlot(AB, idents = 'NSC/Astrocyte', features = c("NFE2L2","ZHX2","VGLL4","ZFP36L2","NR2E1"), cols = c("lightgrey", "red"),split.by = "Age", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

