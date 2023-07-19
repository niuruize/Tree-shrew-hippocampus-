################################################################################
#gene set score
################################################################################
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)

load("~/Downloads/snRNA/cross_species_snRNA.RData")
snRNA$celltype_species <- paste(snRNA$celltype, snRNA$species, sep = "_")
#AD_score
AD_gene <- read_xlsx("AD.xlsx")
gene <- as.list(AD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "AD")
colnames(AB@meta.data)[11]<-"AD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","AD_Score"))
P1<-ggplot(data,aes(celltype.group, AD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "AD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#ADHD_score
ADHD_gene <- read_xlsx("ADHD.xlsx")
gene <- as.list(ADHD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ADHD")
colnames(AB@meta.data)[11]<-"ADHD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","ADHD_Score"))
P1<-ggplot(data,aes(celltype.group, ADHD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ADHD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#ALD_score
ALD_gene <- read_xlsx("ALD.xlsx")
gene <- as.list(ALD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ALD")
colnames(AB@meta.data)[11]<-"ALD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","ALD_Score"))
P1<-ggplot(data,aes(celltype.group, ALD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ALD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#ANO_score
ANO_gene <- read_xlsx("ANO.xlsx")
gene <- as.list(ANO_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ANO")
colnames(AB@meta.data)[11]<-"ANO_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","ANO_Score"))
P1<-ggplot(data,aes(celltype.group, ANO_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ANO_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Anxiety_score
Anxiety_gene <- read_xlsx("Anxiety.xlsx")
gene <- as.list(Anxiety_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Anxiety")
colnames(AB@meta.data)[11]<-"Anxiety_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Anxiety_Score"))
P1<-ggplot(data,aes(celltype.group, Anxiety_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Anxiety_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#ASD_score
ASD_gene <- read_xlsx("ASD.xlsx")
gene <- as.list(ASD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "ASD")
colnames(AB@meta.data)[11]<-"ASD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","ASD_Score"))
P1<-ggplot(data,aes(celltype.group, ASD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "ASD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#BIP_score
BIP_gene <- read_xlsx("BIP.xlsx")
gene <- as.list(BIP_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "BIP")
colnames(AB@meta.data)[11]<-"BIP_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","BIP_Score"))
P1<-ggplot(data,aes(celltype.group, BIP_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "BIP_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Epilepsy_score
Epilepsy_gene <- read_xlsx("Epilepsy.xlsx")
gene <- as.list(Epilepsy_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Epilepsy")
colnames(AB@meta.data)[11]<-"Epilepsy_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Epilepsy_Score"))
P1<-ggplot(data,aes(celltype.group, Epilepsy_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Epilepsy_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Learning_memory_score
Learning_memory_gene <- read_xlsx("Learning_memory.xlsx")
gene <- as.list(Learning_memory_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Learning_memory")
colnames(AB@meta.data)[11]<-"Learning_memory_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Learning_memory_Score"))
P1<-ggplot(data,aes(celltype.group, Learning_memory_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Learning_memory_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#MDD_score
MDD_gene <- read_xlsx("MDD.xlsx")
gene <- as.list(MDD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "MDD")
colnames(AB@meta.data)[11]<-"MDD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","MDD_Score"))
P1<-ggplot(data,aes(celltype.group, MDD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "MDD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Narcolepsy_score
Narcolepsy_gene <- read_xlsx("Narcolepsy.xlsx")
gene <- as.list(Narcolepsy_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Narcolepsy")
colnames(AB@meta.data)[11]<-"Narcolepsy_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Narcolepsy_Score"))
P1<-ggplot(data,aes(celltype.group, Narcolepsy_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Narcolepsy_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#SCZ_score
SCZ_gene <- read_xlsx("SCZ.xlsx")
gene <- as.list(SCZ_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "SCZ")
colnames(AB@meta.data)[11]<-"SCZ_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","SCZ_Score"))
P1<-ggplot(data,aes(celltype.group, SCZ_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "SCZ_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#TOS_score
TOS_gene <- read_xlsx("TOS.xlsx")
gene <- as.list(TOS_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "TOS")
colnames(AB@meta.data)[11]<-"TOS_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","TOS_Score"))
P1<-ggplot(data,aes(celltype.group, TOS_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "TOS_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#OCD_score
OCD_gene <- read_xlsx("OCD.xlsx")
gene <- as.list(OCD_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "OCD")
colnames(AB@meta.data)[11]<-"OCD_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","OCD_Score"))
P1<-ggplot(data,aes(celltype.group, OCD_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "OCD_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Aging_Up_score
Aging_Up_gene <- read_xlsx("Aging_Up.xlsx")
gene <- as.list(Aging_Up_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Aging_Up")
colnames(AB@meta.data)[11]<-"Aging_Up_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Aging_Up_Score"))
P1<-ggplot(data,aes(celltype.group, Aging_Up_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Aging_Up_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Aging_Down_score
Aging_Down_gene <- read_xlsx("Aging_Down.xlsx")
gene <- as.list(Aging_Down_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Aging_Down")
colnames(AB@meta.data)[11]<-"Aging_Down_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Aging_Down_Score"))
P1<-ggplot(data,aes(celltype.group, Aging_Down_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Aging_Down_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

#Aging_HuMi_score
Aging_HuMi_gene <- read_xlsx("Aging_HuMi.xlsx")
gene <- as.list(Aging_HuMi_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Aging_HuMi")
colnames(AB@meta.data)[11]<-"Aging_HuMi_Score" 
##box plot
data<-FetchData(AB, vars = c("celltype_species","Aging_HuMi_Score"))
P1<-ggplot(data,aes(celltype.group, Aging_HuMi_Score))+
  geom_boxplot(aes(fill=celltype.group),outlier.alpha=0,show.legend = F)+theme_bw()+RotatedAxis()
ggsave(filename = "Aging_HuMi_score.pdf", plot = P1, device = 'pdf', width = 50, height = 8, units = 'cm')

