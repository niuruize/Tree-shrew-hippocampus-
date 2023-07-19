#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
#loading data
#Pre-processing and quality control
################################################################################

#CreateSeuratObject---Infancy1
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Infancy1")
Infancy1 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Infancy1", min.cells = 3, min.features = 200)
Infancy1 <- RenameCells(Infancy1, add.cell.id = "Infancy1")
Infancy1[["SampleID"]] <- "Infancy1"
Infancy1[["Age"]] <- "Infancy"
Infancy1[["Sex"]] <- "F"
Infancy1[["percent.mt"]] <- PercentageFeatureSet(Infancy1,pattern = "^MT-")
summary(Infancy1[[]]$percent.mt)
VlnPlot(Infancy1, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Infancy1 <- subset(Infancy1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Infancy1, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Infancy1[[]]$nFeature_RNA)
plotInfancy1_1 <- FeatureScatter(Infancy1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotInfancy1_2 <- FeatureScatter(Infancy1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotInfancy1_1 + plotInfancy1_2
#DoubletFinder
library(DoubletFinder)
Infancy1 <- SCTransform(Infancy1, vars.to.regress = "percent.mt", verbose = FALSE)
Infancy1 <- RunPCA(Infancy1)
Infancy1 <- FindNeighbors(Infancy1, dims = 1:30)
Infancy1 <- FindClusters(Infancy1, resolution = 0.8)
Infancy1 <- RunUMAP(Infancy1, dims = 1:30)
sweep.res.list <- paramSweep_v3(Infancy1, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Infancy1)*8*1e-6
homotypic.prop <- modelHomotypic(Infancy1$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Infancy1)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Infancy1 <- doubletFinder_v3(Infancy1, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Infancy1@meta.data)[ncol(Infancy1@meta.data)]="DoubletFinder"
DimPlot(Infancy1,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Infancy1@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Infancy1 <- CreateSeuratObject(Infancy1@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Infancy2
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Infancy2")
Infancy2 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Infancy2", min.cells = 3, min.features = 200)
Infancy2 <- RenameCells(Infancy2, add.cell.id = "Infancy2")
Infancy2[["SampleID"]] <- "Infancy2"
Infancy2[["Age"]] <- "Infancy"
Infancy2[["Sex"]] <- "M"
Infancy2[["percent.mt"]] <- PercentageFeatureSet(Infancy2,pattern = "^MT-")
summary(Infancy2[[]]$percent.mt)
VlnPlot(Infancy2, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Infancy2 <- subset(Infancy2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Infancy2, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Infancy2[[]]$nFeature_RNA)
plotInfancy2_1 <- FeatureScatter(Infancy2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotInfancy2_2 <- FeatureScatter(Infancy2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotInfancy2_1 + plotInfancy2_2
#DoubletFinder
library(DoubletFinder)
Infancy2 <- SCTransform(Infancy2, vars.to.regress = "percent.mt", verbose = FALSE)
Infancy2 <- RunPCA(Infancy2)
Infancy2 <- FindNeighbors(Infancy2, dims = 1:30)
Infancy2 <- FindClusters(Infancy2, resolution = 0.8)
Infancy2 <- RunUMAP(Infancy2, dims = 1:30)
sweep.res.list <- paramSweep_v3(Infancy2, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Infancy2)*8*1e-6
homotypic.prop <- modelHomotypic(Infancy2$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Infancy2)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Infancy2 <- doubletFinder_v3(Infancy2, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Infancy2@meta.data)[ncol(Infancy2@meta.data)]="DoubletFinder"
DimPlot(Infancy2,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Infancy2@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Infancy2 <- CreateSeuratObject(Infancy2@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Infancy3
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Infancy3")
Infancy3 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Infancy3", min.cells = 3, min.features = 200)
Infancy3 <- RenameCells(Infancy3, add.cell.id = "Infancy3")
Infancy3[["SampleID"]] <- "Infancy3"
Infancy3[["Age"]] <- "Infancy"
Infancy3[["Sex"]] <- "M"
Infancy3[["percent.mt"]] <- PercentageFeatureSet(Infancy3,pattern = "^MT-")
summary(Infancy3[[]]$percent.mt)
VlnPlot(Infancy3, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Infancy3 <- subset(Infancy3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Infancy3, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Infancy3[[]]$nFeature_RNA)
plotInfancy3_1 <- FeatureScatter(Infancy3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotInfancy3_2 <- FeatureScatter(Infancy3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotInfancy3_1 + plotInfancy3_2
#DoubletFinder
library(DoubletFinder)
Infancy3 <- SCTransform(Infancy3, vars.to.regress = "percent.mt", verbose = FALSE)
Infancy3 <- RunPCA(Infancy3)
Infancy3 <- FindNeighbors(Infancy3, dims = 1:30)
Infancy3 <- FindClusters(Infancy3, resolution = 0.8)
Infancy3 <- RunUMAP(Infancy3, dims = 1:30)
sweep.res.list <- paramSweep_v3(Infancy3, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Infancy3)*8*1e-6
homotypic.prop <- modelHomotypic(Infancy3$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Infancy3)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Infancy3 <- doubletFinder_v3(Infancy3, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Infancy3@meta.data)[ncol(Infancy3@meta.data)]="DoubletFinder"
DimPlot(Infancy3,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Infancy3@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Infancy3 <- CreateSeuratObject(Infancy3@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Adult1
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Adult1")
Adult1 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Adult1", min.cells = 3, min.features = 200)
Adult1 <- RenameCells(Adult1, add.cell.id = "Adult1")
Adult1[["SampleID"]] <- "Adult1"
Adult1[["Age"]] <- "Adult"
Adult1[["Sex"]] <- "F"
Adult1[["percent.mt"]] <- PercentageFeatureSet(Adult1,pattern = "^MT-")
summary(Adult1[[]]$percent.mt)
VlnPlot(Adult1, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Adult1 <- subset(Adult1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Adult1, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Adult1[[]]$nFeature_RNA)
plotAdult1_1 <- FeatureScatter(Adult1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotAdult1_2 <- FeatureScatter(Adult1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotAdult1_1 + plotAdult1_2
#DoubletFinder
library(DoubletFinder)
Adult1 <- SCTransform(Adult1, vars.to.regress = "percent.mt", verbose = FALSE)
Adult1 <- RunPCA(Adult1)
Adult1 <- FindNeighbors(Adult1, dims = 1:30)
Adult1 <- FindClusters(Adult1, resolution = 0.8)
Adult1 <- RunUMAP(Adult1, dims = 1:30)
sweep.res.list <- paramSweep_v3(Adult1, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Adult1)*8*1e-6
homotypic.prop <- modelHomotypic(Adult1$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Adult1)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Adult1 <- doubletFinder_v3(Adult1, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Adult1@meta.data)[ncol(Adult1@meta.data)]="DoubletFinder"
DimPlot(Adult1,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Adult1@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Adult1 <- CreateSeuratObject(Adult1@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Adult2
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Adult2")
Adult2 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Adult2", min.cells = 3, min.features = 200)
Adult2 <- RenameCells(Adult2, add.cell.id = "Adult2")
Adult2[["SampleID"]] <- "Adult2"
Adult2[["Age"]] <- "Adult"
Adult2[["Sex"]] <- "M"
Adult2[["percent.mt"]] <- PercentageFeatureSet(Adult2,pattern = "^MT-")
summary(Adult2[[]]$percent.mt)
VlnPlot(Adult2, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Adult2 <- subset(Adult2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Adult2, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Adult2[[]]$nFeature_RNA)
plotAdult2_1 <- FeatureScatter(Adult2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotAdult2_2 <- FeatureScatter(Adult2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotAdult2_1 + plotAdult2_2
#DoubletFinder
library(DoubletFinder)
Adult2 <- SCTransform(Adult2, vars.to.regress = "percent.mt", verbose = FALSE)
Adult2 <- RunPCA(Adult2)
Adult2 <- FindNeighbors(Adult2, dims = 1:30)
Adult2 <- FindClusters(Adult2, resolution = 0.8)
Adult2 <- RunUMAP(Adult2, dims = 1:30)
sweep.res.list <- paramSweep_v3(Adult2, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Adult2)*8*1e-6
homotypic.prop <- modelHomotypic(Adult2$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Adult2)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Adult2 <- doubletFinder_v3(Adult2, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Adult2@meta.data)[ncol(Adult2@meta.data)]="DoubletFinder"
DimPlot(Adult2,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Adult2@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Adult2 <- CreateSeuratObject(Adult2@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Old1
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Old1")
Old1 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Old1", min.cells = 3, min.features = 200)
Old1 <- RenameCells(Old1, add.cell.id = "Old1")
Old1[["SampleID"]] <- "Old1"
Old1[["Age"]] <- "Old"
Old1[["Sex"]] <- "F"
Old1[["percent.mt"]] <- PercentageFeatureSet(Old1,pattern = "^MT-")
summary(Old1[[]]$percent.mt)
VlnPlot(Old1, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Old1 <- subset(Old1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Old1, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Old1[[]]$nFeature_RNA)
plotOld1_1 <- FeatureScatter(Old1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotOld1_2 <- FeatureScatter(Old1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotOld1_1 + plotOld1_2
#DoubletFinder
library(DoubletFinder)
Old1 <- SCTransform(Old1, vars.to.regress = "percent.mt", verbose = FALSE)
Old1 <- RunPCA(Old1)
Old1 <- FindNeighbors(Old1, dims = 1:30)
Old1 <- FindClusters(Old1, resolution = 0.8)
Old1 <- RunUMAP(Old1, dims = 1:30)
sweep.res.list <- paramSweep_v3(Old1, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Old1)*8*1e-6
homotypic.prop <- modelHomotypic(Old1$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Old1)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Old1 <- doubletFinder_v3(Old1, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Old1@meta.data)[ncol(Old1@meta.data)]="DoubletFinder"
DimPlot(Old1,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Old1@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Old1 <- CreateSeuratObject(Old1@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Old2
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Old2")
Old2 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Old2", min.cells = 3, min.features = 200)
Old2 <- RenameCells(Old2, add.cell.id = "Old2")
Old2[["SampleID"]] <- "Old2"
Old2[["Age"]] <- "Old"
Old2[["Sex"]] <- "M"
Old2[["percent.mt"]] <- PercentageFeatureSet(Old2,pattern = "^MT-")
summary(Old2[[]]$percent.mt)
VlnPlot(Old2, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Old2 <- subset(Old2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Old2, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Old2[[]]$nFeature_RNA)
plotOld2_1 <- FeatureScatter(Old2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotOld2_2 <- FeatureScatter(Old2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotOld2_1 + plotOld2_2
#DoubletFinder
library(DoubletFinder)
Old2 <- SCTransform(Old2, vars.to.regress = "percent.mt", verbose = FALSE)
Old2 <- RunPCA(Old2)
Old2 <- FindNeighbors(Old2, dims = 1:30)
Old2 <- FindClusters(Old2, resolution = 0.8)
Old2 <- RunUMAP(Old2, dims = 1:30)
sweep.res.list <- paramSweep_v3(Old2, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Old2)*8*1e-6
homotypic.prop <- modelHomotypic(Old2$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Old2)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Old2 <- doubletFinder_v3(Old2, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Old2@meta.data)[ncol(Old2@meta.data)]="DoubletFinder"
DimPlot(Old2,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Old2@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Old2 <- CreateSeuratObject(Old2@assays$RNA@counts, meta.data = cellinfo)

#CreateSeuratObject---Old3
counts <- Read10X(data.dir = "~/scRNAseq/Hip/rawdata/Old3")
Old3 <- CreateSeuratObject(counts = counts, meta.data = pred.test, project = "Old3", min.cells = 3, min.features = 200)
Old3 <- RenameCells(Old3, add.cell.id = "Old3")
Old3[["SampleID"]] <- "Old3"
Old3[["Age"]] <- "Old"
Old3[["Sex"]] <- "M"
Old3[["percent.mt"]] <- PercentageFeatureSet(Old3,pattern = "^MT-")
summary(Old3[[]]$percent.mt)
VlnPlot(Old3, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Old3 <- subset(Old3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(Old3, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(Old3[[]]$nFeature_RNA)
plotOld3_1 <- FeatureScatter(Old3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotOld3_2 <- FeatureScatter(Old3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotOld3_1 + plotOld3_2
#DoubletFinder
library(DoubletFinder)
Old3 <- SCTransform(Old3, vars.to.regress = "percent.mt", verbose = FALSE)
Old3 <- RunPCA(Old3)
Old3 <- FindNeighbors(Old3, dims = 1:30)
Old3 <- FindClusters(Old3, resolution = 0.8)
Old3 <- RunUMAP(Old3, dims = 1:30)
sweep.res.list <- paramSweep_v3(Old3, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Old3)*8*1e-6
homotypic.prop <- modelHomotypic(Old3$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Old3)))
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Old3 <- doubletFinder_v3(Old3, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Old3@meta.data)[ncol(Old3@meta.data)]="DoubletFinder"
DimPlot(Old3,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
cellinfo <- subset(Old3@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","SampleID","Age","Sex","DoubletFinder"))
Old3 <- CreateSeuratObject(Old3@assays$RNA@counts, meta.data = cellinfo)

################################################################################
#data integration
################################################################################

AB <- merge(Infancy1, y=c(Infancy2,Infancy3,Adult1,Adult2,Old1,Old2,Old3))
AB.list <- SplitObject(AB, split.by = "SampleID")
AB.list <- lapply(X = AB.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000,verbose = FALSE)
})
AB.features <- SelectIntegrationFeatures(object.list = AB.list, nfeatures = 2000)
AB.anchors <- FindIntegrationAnchors(object.list = AB.list, anchor.features = AB.features)
AB <- IntegrateData(anchorset = AB.anchors)
AB
dim(AB[["RNA"]]@counts)
dim(AB[["RNA"]]@data)
dim(AB[["integrated"]]@counts)
dim(AB[["integrated"]]@data)
save(AB,file="scRNA.RData")

rm("Infancy1",'Infancy2',"Infancy3","Adult1","Adult2",'Old1','Old2','Old3'，AB.list,AB.features,AB.anchors)

################################################################################
#Determine the number of clusters in the dataset 
################################################################################

#method 1：Jackstraw
AB <- JackStraw(AB, num.replicate = 100) 
AB <- ScoreJackStraw(AB, dims = 1:20) 
p <- JackStrawPlot(AB, dims = 1:20)
ggsave(filename = "JackStrawPlot.pdf", plot = p6, device = 'pdf', width = 14, height = 12, units = 'cm')
#method 2：Elbow diagram
p <- ElbowPlot(AB)
p
ggsave(filename = "ElbowPlot.pdf", plot = p7, device = 'pdf', width = 14, height = 12, units = 'cm')

################################################################################
#data dimension reduction; find clusters
################################################################################

DefaultAssay(AB) <- "integrated"
# Run the standard workflow for visualization and clustering
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, npcs = 30, verbose = FALSE)
AB <- FindNeighbors(AB, dims = 1:30)
AB <- FindClusters(AB, resolution = 0.8)
AB <- RunUMAP(AB, dims = 1:30)
AB <- RunTSNE(AB, dims = 1:30) 
head(Idents(AB), 5) 
save(AB,file="AB.RData")
