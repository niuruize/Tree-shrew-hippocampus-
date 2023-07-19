
################################################################################
# Supplementary Figure 11
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

load('/Users/niuruize/Downloads/snRNAseq/Hip/10.cross_species/Microglia/cross_species_Micro.RData')

################################################################################
# Supplementary Figure 11B
################################################################################
DefaultAssay(snRNA) <- "RNA"
markers.to.plot <- c("PLD5","CYTL1","BLNK","PIPOX","MEF2C","TOMM40","TIMD4","MAN1A1","FCHSD2","PID1")#top10
pdf("TS_Micro_specific_marker.pdf", width = 15.5,height = 7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "species", dot.scale = 8) +  RotatedAxis()
dev.off()

################################################################################
# Supplementary Figure 11C
################################################################################
markers.to.plot <- c("LRMDA","MAN1A1","ST6GALNAC3","FCHSD2","ELMO1","CAMK1D","HDAC9","GAB2","MAML3") #special
TS_Microglia_specific_marker <- VlnPlot(snRNA,group.by = "group",features = markers.to.plot, stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=90,size=6))
ggsave(filename = "TS_Human_Macaque_Microglia_specific_marker.pdf", plot = Microglia_specific_marker, device = 'pdf', width = 12, height = 20, units = 'cm')
rm('all_marker')

