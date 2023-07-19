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
#human data analysis
################################################################################
counts <-read.table('GSE163737_HHP_gene_expression_matrix.txt', header = TRUE)
Human <- CreateSeuratObject(counts = counts, project = "Human", min.cells = 3, min.features = 200) 
Human[["percent.mt"]] <- PercentageFeatureSet(Human,pattern = "^MT-")
Human <- subset(Human, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 25)
VlnPlot(Human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plotHuman.pdf", width = 28, height = 25, units = "cm")
plotHuman_1 <- FeatureScatter(Human, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotHuman_2 <- FeatureScatter(Human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotHuman_1 + plotHuman_2
ggsave("Human_1.pdf", width = 28, height = 25, units = "cm")
Human.list <- SplitObject(Human, split.by = "orig.ident")
Human.list <- lapply(X = Human.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000,verbose = FALSE)
})
Human.features <- SelectIntegrationFeatures(object.list = Human.list, nfeatures = 2000)
Human.anchors <- FindIntegrationAnchors(object.list = Human.list, anchor.features = Human.features)
Human <- IntegrateData(anchorset = Human.anchors)
DefaultAssay(Human) <- "integrated"
# Run the standard workflow for visualization and clustering
Human <- ScaleData(Human, features = rownames(Human))
Human <- RunPCA(Human, npcs = 30, verbose = FALSE)
Human <- FindNeighbors(Human, dims = 1:30)
Human <- FindClusters(Human, resolution = 0.8)
Human <- RunUMAP(Human, dims = 1:30)
Human <- RunTSNE(Human, dims = 1:30) 
#cell type identification
Idents(Human) <- "seurat_clusters"
Human_all_marker <- VlnPlot(Human,features = c("SLC1A3","ASCL1","HMGB2","SOX2","GFAP","PAX6","HOPX","AQP4","SOX6",
                                       "PTPRC","C1QA","C1QB","CSF1R","VCAN","OLIG2","PDGFRA","BCAS1","FYN",
                                       "MOG","PLP1","ST18","RBFOX3", "SNAP25","PROX1","DCX","MKI67","NEUROD2","SATB2",
                                       "CAMK2A","GAD1","GAD2","LHX6","NR2F2","SLC6A1","CCK","LAMP5","SV2C","CNR1",
                                       "SST","CALB2","PVALB","RELN","EBF1","RGS5"), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6))
ggsave(filename = "Human_all_marker.pdf", plot = Human_all_marker, device = 'pdf', width = 20, height = 20, units = 'cm')
#Specifies the cell type for the cluster
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
                         "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35")
new.cluster.ids <- c("Oligodendrocyte","NSC/Astrocyte","CA1","GABAergic","CA2","GABAergic","NB/ImmN","GABAergic","CA3","Granule","GABAergic", #0-10
                     "Oligodendrocyte","Oligodendrocyte","Microglia","CA3","CA1","Endothelial","Oligodendrocyte","NSC/Astrocyte","Granule","Microglia", #11-20
                     "Oligodendrocyte","Microglia","OPC","OPC","CA2","CA2","GABAergic","Pericyte","NSC/Astrocyte","CA1", #21-30
                     "Ependymal","NB/ImmN","NSC/Astrocyte","Oligodendrocyte","GABAergic") 
Human@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Human@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Human$celltype1 <- factor(Human$celltype1, levels=c("NSC/Astrocyte","NB/ImmN","Granule","CA1","CA2","CA3","GABAergic","OPC","Oligodendrocyte","Microglia","Endothelial","Ependymal","Pericyte"))
table(Human@meta.data$celltype)
save(Human,file="Human.RData")

################################################################################
#Macaque data analysis
################################################################################
Convert('~/GSE163737_RAW/GSE163737_MHP_year4_11_23_raw.h5ad', "h5seurat", overwrite = TRUE,assay = "RNA")
Macaque <- LoadH5Seurat("~/GSE163737_RAW/GSE163737_MHP_year4_11_23_raw.h5seurat")
count = Macaque@assays[["RNA"]]@counts
Macaque <- CreateSeuratObject(counts = count, project = "Macaque", min.cells = 3, min.features = 200) 
MT.genes <- c("COX2","ND5","ATP6","ND3","COX3","COX1","CYTB","ND1","ND6","ND2","ND4","ND4L")
MT.genes <- CaseMatch(MT.genes, rownames(Macaque))
Macaque[["percent.mt"]] <- PercentageFeatureSet(Macaque,features = MT.genes)
summary(Macaque[[]]$percent.mt)
Macaque <- subset(Macaque, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
VlnPlot(Macaque, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plotMacaque.pdf", width = 28, height = 25, units = "cm")
plotMacaque_1 <- FeatureScatter(Macaque, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotMacaque_2 <- FeatureScatter(Macaque, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMacaque_1 + plotMacaque_2
ggsave("Macaque_1.pdf", width = 28, height = 25, units = "cm")
Macaque <- NormalizeData(Macaque, normalization.method = "LogNormalize", scale.factor = 10000)
Macaque <- NormalizeData(Macaque)
Macaque <- FindVariableFeatures(Macaque, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Macaque), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Macaque)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(Macaque)
Macaque <- ScaleData(Macaque, features = all.genes)
Macaque <- RunPCA(Macaque, features = VariableFeatures(object = Macaque))
Macaque <- JackStraw(Macaque, num.replicate = 100)
Macaque <- ScoreJackStraw(Macaque, dims = 1:50)
JackStrawPlot(Macaque, dims = 1:30)
ElbowPlot(Macaque)
Macaque <- FindNeighbors(Macaque, dims = 1:30)
Macaque <- FindClusters(Macaque, resolution = 0.8)
Macaque <- RunUMAP(Macaque, dims = 1:30)
Macaque = Macaque[,Macaque@meta.data[["seurat_clusters"]] %in% c("0","1","2","3","5","6","7","8","9","10","11","12","13","14","15","16","18","19","20",
                         "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")]
current.cluster.ids <- c("0","1","2","3","5","6","7","8","9","10","11","12","13","14","15","16","18","19","20",
                         "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
new.cluster.ids <- c("NFOL","Unk","DG_ExN","nonDG_ExN","nonDG_ExN","OPC","InN","nonDG_ExN","NSC","InN", #0-10
                     "InN","Astroglia","Microglia","nonDG_ExN","ImN","InN","nonDG_ExN","nonDG_ExN","nonDG_ExN", #11-20
                     "MOL","nonDG_ExN","OPC","nonDG_ExN","nonDG_ExN","nonDG_ExN","Microglia","InN","ImN","nonDG_ExN", #21-30
                     "InN","NFOL","Astroglia","nonDG_ExN","NSC","OPC","OPC","MOL","End","NSC") 
Macaque@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Macaque@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Macaque@meta.data$celltype)
Macaque@meta.data[["celltype"]]<-factor(Macaque@meta.data[["celltype"]], levels=c("NSC","Astroglia","Microglia","OPC","NFOL","ImN","DG_ExN","nonDG_ExN","InN","End","Unk"))
#Macaque to Human
library(biomaRt)
listMarts()
umart = useMart('ensembl')
datalist = listDatasets(umart)
head(datalist)
searchDatasets(mart = umart, pattern = "Tree")
macaque <- useMart("ensembl", dataset = "mfascicularis_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
Human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
Macaque_Human <- getLDS(attributes = c("external_gene_name"),filters = "external_gene_name",values = rownames(Macaque),mart = macaque,
                   attributesL = c("external_gene_name"),martL = Human,uniqueRows = T)
write.table(Macaque_Human,file="Macaque_Human.txt",quote=F,sep="\t",row.names=F,col.names=T)
Macaque <- subset(Macaque, features = Macaque_Human$Gene.name)
rownames(Macaque) <- Macaque_Human$Gene.name.1
save(Macaque,file="Macaque.RData")

################################################################################
#Mouse data analysis
################################################################################
counts <- read.table('~/GSE104323_10X_expression_data_V2.tab.gz', row.names = 1,header = T, sep = "\t")
barcodes <- read.table('~/GSE104323_metadata_barcodes_24185cells.txt', header = T, sep = "\t")
rownames(barcodes) <- make.unique(barcodes$Sample.name..24185.single.cells.)
barcodes <- data.frame(barcodes[,2:11])
barcodes <- data.frame(barcodes[1:24185,])
colnames(counts) <- rownames(barcodes)
Mouse <- CreateSeuratObject(counts = counts, assay = 'RNA',meta.data = barcodes, project = "Mouse", min.cells = 3, min.features = 200)
Mouse[["percent.mt"]] <- PercentageFeatureSet(Mouse,pattern = "^mt-")
Mouse <- subset(Mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 25)
VlnPlot(Mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plotMouse.pdf", width = 28, height = 25, units = "cm")
plotMouse_1 <- FeatureScatter(Mouse, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotMouse_2 <- FeatureScatter(Mouse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMouse_1 + plotMouse_2
ggsave("Mouse_1.pdf", width = 28, height = 25, units = "cm")
Mouse.list <- SplitObject(Mouse, split.by = "orig.ident")
Mouse.list <- lapply(X = Mouse.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000,verbose = FALSE)
})
Mouse.features <- SelectIntegrationFeatures(object.list = Mouse.list, nfeatures = 2000)
Mouse.anchors <- FindIntegrationAnchors(object.list = Mouse.list, anchor.features = Mouse.features)
Mouse <- IntegrateData(anchorset = Mouse.anchors)
DefaultAssay(Mouse) <- "integrated"
# Run the standard workflow for visualization and clustering
Mouse <- ScaleData(Mouse, features = rownames(Mouse))
Mouse <- RunPCA(Mouse, npcs = 30, verbose = FALSE)
Mouse <- FindNeighbors(Mouse, dims = 1:30)
Mouse <- FindClusters(Mouse, resolution = 0.8)
Mouse <- RunUMAP(Mouse, dims = 1:30)
Mouse <- RunTSNE(Mouse, dims = 1:30) 
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10",
                         "11","12","13","14","15","16","17","18","19","20",
                         "21","22","23","24","25","26","27","28","29","30",
                         "31","32","33","34","35","36","37","38","39","42")
new.cluster.ids <- c("MOL","GC","GC","Pyr","OPC","Pyr","GC","Astrocyte","MOL","GC_im","NB", #0-10
                     "NB","Pyr","GABA_im","Astrocyte","Pyr_im","Pyr_im","Microglia","RGL","GABA","GC",#11-20
                     "GABA_im", "GABA_im","nIPC_p","Pyr","Endothelial","nIPC","GABA_im","Pyr_im","NFOL","GC",#21-29
                     "Pyr_im","GABA","Ependymal","GC","VLMC","RGL","CR","RGL","Pyr","Pyr") 
Mouse@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Mouse@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Mouse@meta.data[["celltype"]]<-factor(Mouse@meta.data[["celltype"]], levels=c("RGL","Astrocyte","nIPC_p","nIPC","NB","GC_im","GC","Pyr_im","Pyr","GABA_im","GABA","OPC","NFOL","MOL","Microglia","CR","Ependymal","Endothelial","VLMC"))
table(Mouse@meta.data$celltype)
#Mouse to Human
library(biomaRt)
listMarts()
umart = useMart('ensembl')
datalist = listDatasets(umart)
head(datalist)
searchDatasets(mart = umart, pattern = "Tree")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
Human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
Mouse_Human <- getLDS(attributes = c("external_gene_name"),filters = "external_gene_name",values = rownames(Mouse),mart = mouse,
                     attributesL = c("external_gene_name"),martL = Human,uniqueRows = T)
write.table(Mouse_Human,file="Mouse_Human.txt",quote=F,sep="\t",row.names=F,col.names=T)
Mouse <- subset(Mouse, features = Mouse_Human$Gene.name)
rownames(Mouse) <- Mouse_Human$Gene.name.1
save(Mouse,file="Mouse.RData")

################################################################################
#TS data analysis
################################################################################
load("~/TS.RData")
#TS to Human
library(biomaRt)
listMarts()
umart = useMart('ensembl')
datalist = listDatasets(umart)
head(datalist)
searchDatasets(mart = umart, pattern = "Tree")
TS <- useMart('ensembl',dataset = "tbelangeri_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
Human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
TS_Human <- getLDS(attributes = c("external_gene_name"),filters = "external_gene_name",values = rownames(TS),mart = TS,
                     attributesL = c("external_gene_name"),martL = Human,uniqueRows = T)
write.table(TS_Human,file="TS_Human.txt",quote=F,sep="\t",row.names=F,col.names=T)
TS <- subset(TS, features = TS_Human$Gene.name)
rownames(TS) <- TS_Human$Gene.name.1
save(TS,file="TS.RData")

################################################################################
#cross_species_integration
################################################################################
#loading data
load("~/TS.RData")
load("~/Human.RData")
load("~/Mouse.RData")
load("~/Macaque.RData")

TS@meta.data$species="TS"
Human@meta.data$species="Human"
Mouse@meta.data$species="Mouse"
Macaque@meta.data$species="Macaque"

################################################################################
#data integration
################################################################################
snRNA <- merge(TS, y=c(Human, Mouse, Macaque))
table(snRNA$orig.ident)
cellinfo <- subset(snRNA@meta.data, select= c("orig.ident","percent.mt","celltype","Age","species"))
snRNA <- CreateSeuratObject(snRNA@assays$RNA@counts, meta.data = cellinfo)
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
snRNA <- IntegrateData(snRNA.anchors, normalization.method = "LogNormalize", dims = 1:50)
snRNA
dim(snRNA[["RNA"]]@counts)
dim(snRNA[["RNA"]]@data)
dim(snRNA[["integrated"]]@counts)
dim(snRNA[["integrated"]]@data)

################################################################################
#data dimension reduction, UMAP and TSNE
################################################################################
DefaultAssay(snRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
snRNA <- ScaleData(snRNA, features = rownames(snRNA))
snRNA <- RunPCA(snRNA, npcs = 100, verbose = FALSE)
snRNA <- FindNeighbors(snRNA, dims = 1:20)
snRNA <- FindClusters(snRNA, resolution = 0.8)
snRNA <- RunUMAP(snRNA, dims = 1:20)
snRNA <- RunTSNE(snRNA, dims = 1:20)
save(snRNA,file="cross_species_snRNA.RData")



