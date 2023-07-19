
################################################################################
# Supplementary Figure 12
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
# Supplementary Figure 12A
################################################################################
load("~/Downloads/snRNA/TS_Micro.RData")
#HuMi_score
HuMi_gene <- read_xlsx("HuMi.xlsx")
gene <- as.list(HuMi_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "HuMi")
colnames(AB@meta.data)[11]<-"HuMi_Score" 
#boxplot
b<-FetchData(AB, vars = c("celltype","HuMi_Score"))
my_comparisons <- list( c("Micro1", "Micro2"), c("Micro1", "Micro3"), c("Micro2", "Micro3"))
ggboxplot(b, x = "celltype", y = "HuMi_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "celltype", palette = "npg")+stat_compare_means(comparisons = my_comparisons)
ggsave("HuMi_subtype.pdf",width = 6,height = 8,units = "cm")

################################################################################
# Supplementary Figure 12 B
################################################################################
my_comparisons <- list( c("Infancy", "Adult"), c("Adult", "Old"), c("Infancy", "Old"))
b<-FetchData(AB, vars = c("orig.ident","HuMi_Score","celltype"))
ggboxplot(b, x = "orig.ident", y = "HuMi_Score",facet.by = "celltype",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "orig.ident", palette = "npg",ncol=3)+stat_compare_means(comparisons = my_comparisons)
ggsave("HuMi_subtype_age.pdf",width = 10,height = 8,units = "cm")

################################################################################
# Supplementary Figure 12 D
################################################################################
# Wnt_score
Wnt_gene <- read_xlsx("Wnt.xlsx")
gene <- as.list(Wnt_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Wnt")
colnames(AB@meta.data)[11]<-"Wnt_Score" 
b<-FetchData(AB, vars = c("orig.ident","Wnt_Score"))
b[["orig.ident"]]<-factor(b[["orig.ident"]], levels=c("Infancy","Adult","Old"))
P1<-ggdensity(b, x = "Wnt_Score",
          add = "mean", rug = TRUE,
          color = "orig.ident", fill = "orig.ident",
          palette = c("#00AFBB", "#E7B800","#8B658B"))
ggsave(filename = "Wnt.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')

################################################################################
# Supplementary Figure 12 E
################################################################################
load("~/Downloads/snRNA/TS_Micro.RData")
markers.to.plot <- c("SALL1","P2RY12","MEF2C","TOMM40","CD74","TLR4","NFKB1","C1QB","P2RY6","TYROBP","IL2RA","IL18","IL4R","IL6ST","TGFBR1","ILRUN","ITGAM","NRROS","SLC7A7")
pdf("Inflammatory_1.pdf", width = 8,height = 2)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) + RotatedAxis()
dev.off()
pdf("Inflammatory_2.pdf", width = 8,height = 2)
DotPlot(snRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "age", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

load("~/Downloads/snRNA/TS_Micro1.RData")
markers.to.plot <- c("SALL1","P2RY12","MEF2C","TOMM40","CD74","TLR4","NFKB1","C1QB","P2RY6","TYROBP","IL2RA","IL18","IL4R","IL6ST","TGFBR1","ILRUN","ITGAM","NRROS","SLC7A7")
pdf("Inflammatory_3.pdf", width = 8,height = 2)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "subcelltype", dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# Supplementary Figure 12F-J
################################################################################

# The cluster-specific genes obtained from Figure 7L were used for gene sets scoring

