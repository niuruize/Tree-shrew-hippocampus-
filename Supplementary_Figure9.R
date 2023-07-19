################################################################################
#---Supplementary Figure 9
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

################################################################################
#---Supplementary Figure 9A
################################################################################
load("~/Downloads/snRNA/cross_species_NSC.RData")
markers.to.plot <- c("SLC1A3","GFAP","AQP4","HES5","SOX2","EMX2","FOXG1","EGFR","LPAR1","ETNPPL",
                     "MFGE8","SOX5","SOX6","CDK1","ASCL1","TFAP2C","EOMES","IGFBPL1","CALB2","PLK5",
                     "NES","NEUROG2","MCM2","DCX","NEUROD1","CALB1","TOP2A","ELAVL2","VIM","NOTCH2",
                     "SOX4","SOX11", "PADI2","FRZB","WNT8B")
pdf("NSC_marker.pdf", width = 15.5,height = 7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) +  RotatedAxis()
dev.off()

################################################################################
#---Supplementary Figure 9B
################################################################################
#RGL1
snRNA1<- snRNA[,snRNA@meta.data$celltype %in% c("RGL1")]
av<-AverageExpression(snRNA1,group.by = "species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_RGL1_group.csv") #保存结果
pdf("cor_RGL1_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果
#RGL2
snRNA2<- snRNA[,snRNA@meta.data$celltype %in% c("RGL2")]
av<-AverageExpression(snRNA2,group.by = "species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_RGL2.csv") #保存结果
pdf("cor_RGL2_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果
#Astrocyte
snRNA2<- snRNA[,snRNA@meta.data$celltype %in% c("Astrocyte")]
av<-AverageExpression(snRNA2,group.by = "species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_Astrocyte.csv") #保存结果
pdf("cor_Astrocyte_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果
#nIPC_p
snRNA2<- snRNA[,snRNA@meta.data$celltype %in% c("nIPC_p")]
av<-AverageExpression(snRNA2,group.by = "NSC.group", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_nIPC_p.csv") #保存结果
pdf("cor_nIPC_p_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果
#NB
snRNA2<- snRNA[,snRNA@meta.data$celltype %in% c("NB")]
av<-AverageExpression(snRNA2,group.by = "species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_NB.csv") #保存结果
pdf("cor_NB_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果
#Immu_RGL
snRNA2<- snRNA[,snRNA@meta.data$celltype %in% c("Immu_RGL")]
av<-AverageExpression(snRNA2,group.by = "species", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_Immu_RGL.csv") #保存结果
pdf("cor_Immu_RGL_group_1.pdf", width = 4.5,height = 4)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_group.csv") #保存结果

################################################################################
#---Supplementary Figure 9C
################################################################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db")

#GO enrichment analysis was performed for each cell type top100 genes separately
rt=read.table("RGL1.txt",sep="\t",check.names=F,header=F)
rt=read.table("RGL2.txt",sep="\t",check.names=F,header=F)
rt=read.table("NB.txt",sep="\t",check.names=F,header=F)
rt=read.table("Astrocyte.txt",sep="\t",check.names=F,header=F)
rt=read.table("nIPC.txt",sep="\t",check.names=F,header=F)
rt=read.table("Immu_RGL.txt",sep="\t",check.names=F,header=F)
rt=read.table("TS_RGL.txt",sep="\t",check.names=F,header=F)

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
id=cbind(rt,1,entrezID=entrezIDs)
colnames(id)=c("symbol","logFC","entrezID")

rt=id
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichGO(gene = gene,ont = "BP",OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, readable =T)
write.csv(kk@result,"RGL1_BP.csv", row.names = F)

#plot the heatmap using TBtools

################################################################################
#---Supplementary Figure 9D
################################################################################
# DEGs between different species within the same cell type:
{snRNA$species_subcelltype <- paste(snRNA$species,snRNA$subcelltype, sep = "_")
table(snRNA$species_subcelltype)
Idents(snRNA)="species_subcelltype"
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Human_RGL1"), ident.2 = ("Mouse_RGL1"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL1_Human_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Macaque_RGL1"), ident.2 = ("Mouse_RGL1"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL1_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("TS_RGL1"), ident.2 = ("Mouse_RGL1"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL1_TS_Mouse.csv")

age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Human_RGL2"), ident.2 = ("Mouse_RGL2"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL2_Human_Mouse.csv")
age.response <- FindMarkers(snRNA,logfc.threshold = 0.01, ident.1 = c("Macaque_RGL2"), ident.2 = ("Mouse_RGL2"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL2_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA,logfc.threshold = 0.01, ident.1 = c("TS_RGL2"), ident.2 = ("Mouse_RGL2"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="RGL2_TS_Mouse.csv")

age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Human_Astrocyte"), ident.2 = ("Mouse_Astrocyte"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Astrocyte_Human_Mouse.csv")
age.response <- FindMarkers(snRNA,logfc.threshold = 0.01, ident.1 = c("Macaque_Astrocyte"), ident.2 = ("Mouse_Astrocyte"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Astrocyte_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("TS_Astrocyte"), ident.2 = ("Mouse_Astrocyte"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Astrocyte_TS_Mouse.csv")

age.response <- FindMarkers(snRNA,logfc.threshold = 0.01, ident.1 = c("Human_nIPC_p"), ident.2 = ("Mouse_nIPC_p"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="nIPC_p_Human_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Macaque_nIPC_p"), ident.2 = ("Mouse_nIPC_p"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="nIPC_p_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("TS_nIPC_p"), ident.2 = ("Mouse_nIPC_p"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="nIPC_p_TS_Mouse.csv")

age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Human_NB"), ident.2 = ("Mouse_NB"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="NB_Human_Mouse.csv")
age.response <- FindMarkers(snRNA,logfc.threshold = 0.01, ident.1 = c("Macaque_NB"), ident.2 = ("Mouse_NB"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="NB_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("TS_NB"), ident.2 = ("Mouse_NB"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="NB_TS_Mouse.csv")

age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Human_Immu_RGL"), ident.2 = ("Mouse_Immu_RGL"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Immu_RGL_Human_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("Macaque_Immu_RGL"), ident.2 = ("Mouse_Immu_RGL"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Immu_RGL_Macaque_Mouse.csv")
age.response <- FindMarkers(snRNA, logfc.threshold = 0.01,ident.1 = c("TS_Immu_RGL"), ident.2 = ("Mouse_Immu_RGL"),test.use = "MAST", verbose = TRUE)
head(age.response, n = 10)
write.csv(age.response,file="Immu_RGL_TS_Mouse.csv")
}
#GO enrichement analysis for human, macaque and TS specific-genes for NB and astrocyte
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db")

rt=read.table("Astrocyte_specific_genes.txt",sep="\t",check.names=F,header=F)
rt=read.table("NB_specific_genes.txt",sep="\t",check.names=F,header=F)

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
id=cbind(rt,1,entrezID=entrezIDs)
colnames(id)=c("symbol","logFC","entrezID")

rt=id
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichGO(gene = gene,ont = "BP",OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, readable =T)
write.csv(kk@result,"NB_BP.csv", row.names = F)

################################################################################
#---Supplementary Figure 9E
################################################################################
markers.to.plot <- c("KCNIP4","ROBO2","STXBP5L","CCSER1","FGF14","NKAIN2","PPFIA1","RGS7","TRPM3","JAZF1","SH3D19", # Human, macaque and TS specific-genes for NB
                     "SAMD4A","GRIP1","DPP10","ZNF704","ENOX1","TOX","PIP4K2A","CCSER1","CACNA2D1","TRPM3","LRRC4C",# Human, macaque and TS specific-genes for Astrocyte
                     "FLRT2", "RGS7","FHIT","AFG1L","CA8","TMEM178B","COL5A3","DLG2","NCAM2","UST","PEX5L","PPFIA2",# Human, macaque and TS specific-genes for Astrocyte
                     "RUNX1T1","DGKB","MOXD1","PRR16","SH3D19","TPD52L1","ARL15","XKR6","TMEM117","MTHFD1L","STXBP5L","RETREG1","SAMD12",) # Human, macaque and TS specific-genes for Astrocyte
pdf("Human_macaque_TS_Astro_NB_marker.pdf", width = 15.5,height = 7)
DotPlot(snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) +  RotatedAxis()
dev.off()


