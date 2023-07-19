################################################################################
# Supplementary Figure 4
# gene set score
################################################################################
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)

################################################################################
# Supplementary Figure 4A,C
################################################################################
load("~/Downloads/snRNA/snRNA.RData")

#AD_score
AD_gene <- read_xlsx("AD.xlsx")
gene <- as.list(AD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "AD")
colnames(AB@meta.data)[11]<-"AD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","AD_Score"))
P1<-ggplot(data,aes(celltype.group, AD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "AD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","AD_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "AD_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("AD_aging.pdf",width = 8,height = 10,units = "cm")

#ADHD_score
ADHD_gene <- read_xlsx("ADHD.xlsx")
gene <- as.list(ADHD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ADHD")
colnames(AB@meta.data)[11]<-"ADHD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","ADHD_Score"))
P1<-ggplot(data,aes(celltype.group, ADHD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ADHD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","ADHD_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "ADHD_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("ADHD_aging.pdf",width = 8,height = 10,units = "cm")

#ANO_score
ANO_gene <- read_xlsx("ANO.xlsx")
gene <- as.list(ANO_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ANO")
colnames(AB@meta.data)[11]<-"ANO_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","ANO_Score"))
P1<-ggplot(data,aes(celltype.group, ANO_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ANO_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","ANO_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "ANO_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("ANO_aging.pdf",width = 8,height = 10,units = "cm")

#ASD_score
ASD_gene <- read_xlsx("ASD.xlsx")
gene <- as.list(ASD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ASD")
colnames(AB@meta.data)[11]<-"ASD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","ASD_Score"))
P1<-ggplot(data,aes(celltype.group, ASD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ASD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","ASD_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "ASD_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("ASD_aging.pdf",width = 8,height = 10,units = "cm")

#BIP_score
BIP_gene <- read_xlsx("BIP.xlsx")
gene <- as.list(BIP_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "BIP")
colnames(AB@meta.data)[11]<-"BIP_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","BIP_Score"))
P1<-ggplot(data,aes(celltype.group, BIP_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "BIP_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","BIP_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "BIP_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("BIP_aging.pdf",width = 8,height = 10,units = "cm")

#MDD_score
MDD_gene <- read_xlsx("MDD.xlsx")
gene <- as.list(MDD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "MDD")
colnames(AB@meta.data)[11]<-"MDD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","MDD_Score"))
P1<-ggplot(data,aes(celltype.group, MDD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "MDD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","MDD_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "MDD_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("MDD_aging.pdf",width = 8,height = 10,units = "cm")

#SCZ_score
SCZ_gene <- read_xlsx("SCZ.xlsx")
gene <- as.list(SCZ_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "SCZ")
colnames(AB@meta.data)[11]<-"SCZ_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","SCZ_Score"))
P1<-ggplot(data,aes(celltype.group, SCZ_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "SCZ_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","SCZ_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "SCZ_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("SCZ_aging.pdf",width = 8,height = 10,units = "cm")

#TOS_score
TOS_gene <- read_xlsx("TOS.xlsx")
gene <- as.list(TOS_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "TOS")
colnames(AB@meta.data)[11]<-"TOS_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","TOS_Score"))
P1<-ggplot(data,aes(celltype.group, TOS_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "TOS_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","TOS_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "TOS_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("TOS_aging.pdf",width = 8,height = 10,units = "cm")

#OCD_score
OCD_gene <- read_xlsx("OCD.xlsx")
gene <- as.list(OCD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "OCD")
colnames(AB@meta.data)[11]<-"OCD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype","OCD_Score"))
P1<-ggplot(data,aes(celltype.group, OCD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "OCD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')
# violin
b<-FetchData(AB, vars = c("orig.ident","OCD_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "OCD_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("OCD_aging.pdf",width = 8,height = 10,units = "cm")

################################################################################
# Supplementary Figure 4B
################################################################################
library(readxl)
genes <- read_xlsx("score_DEGs.xlsx")
p <- VlnPlot(AB, features = c(genes$gene), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=4))
ggsave(filename = "score_DEGs_1.pdf", plot = p, device = 'pdf', width = 10, height = 20, units = 'cm')







