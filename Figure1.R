#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
###---Figure 1B
################################################################################
load("~/snRNA.RData")
p <- DimPlot(snRNA, reduction = "umap", group.by = "celltype", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
ggsave(filename = "celltype_umap.pdf", plot = p, device = 'pdf', width = 10, height = 8, units = 'cm')

################################################################################
###---Figure 1C
################################################################################
# The input data were obtained from cross "species integration.R"
load("~/cross_species_snRNA.RData")
p <- DimPlot(snRNA, reduction = "umap", group.by = "celltype", split.by = "species", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "celltype_cross_species_umap.pdf", plot = p, device = 'pdf', width = 22.5, height =11, units = 'cm')

################################################################################
###---Figure 1D
################################################################################
load("~/cross_species_snRNA.RData")
snRNA$celltype_species <- paste(snRNA$celltype, snRNA$species, sep = "_")
av<-AverageExpression(snRNA,group.by = "celltype_species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype_species.csv") #保存结果
pdf("cor_celltype_species.pdf", width = 20,height = 20)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()












