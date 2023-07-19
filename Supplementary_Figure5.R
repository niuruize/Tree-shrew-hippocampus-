################################################################################
# Supplementary Figure 5A
# Analysis of gene expression trends during hippocampal aging
################################################################################

library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)

load("~/Downloads/snRNA/snRNA.RData")
#bulk
age.averages <- AverageExpression(snRNA,group.by = "Age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["integrated"]]
mat <- as.matrix(df1)
head(mat)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
head(df.s)
m1 <- mestimate(dt.s)
set.seed(007)
cl <- mfuzz(dt.s,c=6,m=m1)

################################################################################
# line chart
################################################################################
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s))
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("SteelBlue","yellow","red")
mycolor <- colorRampPalette(mycol)(50)
cl$size
head(cl$cluster)
head(cl$membership)
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
head(gene_cluster)
#output
write.csv(gene_cluster, 'gene_cluster.csv')

################################################################################
# heatmap
################################################################################
library(pheatmap)
age.averages <- AverageExpression(AB,group.by = "age")
df1 <- data.age.averages[["integrated"]]
mat <- as.matrix(df1)
library(readxl)
markers <- read_xlsx("gene_cluster.xlsx")
mat1 = mat[markers$gene,]
pdf("expression_heatmap_1.pdf", width = 8,height = 6)
pheatmap(mat1,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,scale = "row",show_rownames =F,
         annotation_col = annotation_col,annotation_colors = ann_colors,annotation_names_col = F)
dev.off()

################################################################################
# AddModuleScore
################################################################################
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)

cluster1_gene <- read_xlsx("cluster1_gene.xlsx")
gene<-as.list(cluster1_gene)
AB<-AddModuleScore(snRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[14]<-"cluster1_Score"
#
P1<-VlnPlot(AB, features = 'cluster1_Score')
ggsave(filename = "cluster1_1.pdf", plot = P1, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('P1') 
#
b<-FetchData(AB, vars = c("orig.ident","cluster1_Score"))
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggviolin(b, x = "orig.ident", y = "cluster1_Score", fill = "orig.ident",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("cluster1_2.pdf",width = 8,height = 8,units = "cm")

################################################################################
# GO
################################################################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db")

rt=read.table("cluster1.txt",sep="\t",check.names=F,header=F)

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
id=cbind(rt,1,entrezID=entrezIDs)
colnames(id)=c("symbol","logFC","entrezID")
rt=id
rt=subset(rt,rt$entrezID != "NA")
gene = rt$entrezID
kk <- enrichGO(gene = gene,ont = "BP",OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, readable =T)
write.csv(kk@result,"cluater1.csv", row.names = F)

################################################################################
# Supplementary Figure 5B,C
# hdWGCNA Analysis
################################################################################
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

theme_set(theme_cowplot())
set.seed(12345)
DimPlot(snRNA, group.by = 'celltype', label = T) + umap_theme()

# select genes: fraction >= 5%
snRNA <- SetupForWGCNA(snRNA, gene_select = "fraction", fraction = 0.05, wgcna_name = "tutorial")
length(snRNA@misc$tutorial$wgcna_genes)
# metacells
snRNA <- MetacellsByGroups(snRNA, group.by = c("celltype", "Age"), k = 30, max_shared = 10, ident.group = 'celltype')
snRNA <- NormalizeMetacells(snRNA)
metacell_obj <- GetMetacellObject(snRNA)
snRNA <- NormalizeMetacells(snRNA)
snRNA <- ScaleMetacells(snRNA, features = rownames(snRNA@assays$integrated))
snRNA <- RunPCAMetacells(snRNA, features = rownames(snRNA@assays$integrated))
snRNA <- RunUMAPMetacells(snRNA, dim=1:20)

p1 <- DimPlotMetacells(snRNA, group.by = 'celltype', label = T) + umap_theme()+ggtitle("celltype")
p2 <- DimPlotMetacells(snRNA, group.by = 'Age', label = T) + umap_theme()+ggtitle("celltype")
p1|p2

# Co-expression Network analysis
snRNA <- SetDatExpr(snRNA, group_name = c("Infancy","Adult","Old"), group.by = 'Age', assay = 'RNA', slot = 'data')
snRNA <- TestSoftPowers(snRNA, networkType = 'signed')
plot_list <- PlotSoftPowers(snRNA)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(snRNA)
head(power_table)
# Build co-expression network
snRNA <- ConstructNetwork(snRNA, soft_power = 5, setDatExpr = F, c("Infancy","Adult","Old"))
PlotDendrogram(snRNA, main = "InN hdWGCNA Dendrogram")
snRNA@misc$tutorial$wgcna_modules %>% head
table(snRNA@misc$tutorial$wgcna_modules$module)
TOM <- GetTOM(snRNA)

# Calculate module feature genes
# need to run ScaleData first or else harmony throws an error:
snRNA <- ScaleData(snRNA, features=rownames(snRNA@assays$integrated))

# compute all MEs in the full single-cell dataset
snRNA <- ModuleEigengenes(snRNA, group.by.vars="Age")

# harmonized module eigengenes:
hMEs <- GetMEs(snRNA)
head(hMEs)

# module eigengenes:
MEs <- GetMEs(snRNA, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
snRNA <- ModuleConnectivity(snRNA, group.by = 'age', group_name = c("Infancy","Adult","Old"))

# rename the modules
snRNA <- ResetModuleNames(snRNA, new_name = "M")

# plot genes ranked by kME for each module
pdf("ranked by kME_1.pdf", width = 10,height = 4)
PlotKMEs(snRNA, ncol=4)
dev.off()

# get the module assignment table:
modules <- GetModules(snRNA)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(snRNA, n_hubs = 100)
head(hub_df)
write.csv(hub_df,"hub_df.csv",row.names = F)
saveRDS(snRNA, file='hdWGCNA_object.rds')

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
snRNA <- ModuleExprScore(snRNA,n_genes = 25, method='Seurat')

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
#library(remotes)
#remotes::install_github("carmonalab/UCell")
#library(UCell)
snRNA <- ModuleExprScore(snRNA, n_genes = 25, method='Seurat')

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(snRNA,features='hMEs', # plot the hMEs
                               order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf("cluster_marker_1.pdf", width = 10,height = 6)
wrap_plots(plot_list, ncol=4)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(snRNA, features='scores', # plot the hub gene scores
                               order='shuffle', # order so cells are shuffled
                               ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf("UMAP_2.pdf", width = 10,height = 6)
wrap_plots(plot_list, ncol=4)
dev.off()

# plot module correlagram
ModuleCorrelogram(snRNA)

# get hMEs from seurat object
MEs <- GetMEs(snRNA, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
snRNA@meta.data <- cbind(snRNA@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(snRNA, features=mods, group.by = 'Age') + coord_flip() + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue')
ggsave(filename = "Age_dot.pdf", plot = p, device = 'pdf', width = 10, height = 7, units = 'cm')
p <- DotPlot(snRNA, features=mods, group.by = 'celltype') + coord_flip() + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue')
ggsave(filename = "celltype_dot.pdf", plot = p, device = 'pdf', width = 18, height = 8, units = 'cm')

par(mfrow = c(8,1))
p1 <- VlnPlot(snRNA, features = 'M1',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p2 <- VlnPlot(snRNA, features = 'M2',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p3 <- VlnPlot(snRNA, features = 'M3',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p4 <- VlnPlot(snRNA, features = 'M4',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p5 <- VlnPlot(snRNA, features = 'M5',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p6 <- VlnPlot(snRNA, features = 'M6',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p7 <- VlnPlot(snRNA, features = 'M7',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
p8 <- VlnPlot(snRNA, features = 'M8',group.by = 'celltype',pt.size = 0)+ geom_boxplot(width=.25, fill='white')+ xlab('') + ylab('hME') + NoLegend()
pdf("hME_celltype.pdf", width = 10,height = 20)
p1+p2+p3+p4+p5+p6+p7+p8
dev.off()

# Network Visualization
library(igraph)
ModuleNetworkPlot(snRNA)

# hubgene network
HubGeneNetworkPlot(snRNA, n_hubs = 5, n_other=10, edge_prop = 0.75, mods = 'all')
HubGeneNetworkPlot(snRNA,  return_graph=F)

snRNA <- RunModuleUMAP(snRNA,
                       n_hubs = 10, # number of hub genes to include for the UMAP embedding
                       n_neighbors=15, # neighbors parameter for UMAP
                       min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(snRNA)
# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
ModuleUMAPPlot(snRNA,
               edge.alpha=0.25,
               sample_edges=TRUE,
               edge_prop=0.1, # proportion of edges to sample (20% here)
               label_hubs=5 ,# how many hub genes to plot per module?
               keep_grey_edges=FALSE)
g <- ModuleUMAPPlot(snRNA,  return_graph=TRUE)

#
library(stringr)
MEs <- GetMEs(snRNA, harmonized=FALSE)
MEs_col = MEs
colnames(MEs_col)= paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs), "ME", ""))))
MEs_col = orderMEs(MEs_col)
head(MEs_col)
pdf("hMEs_cor.pdf", width = 4,height = 6.5)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmao", marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), plotDendrograms =T, xLabelsAngle =90)
dev.off()






