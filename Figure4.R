################################################################################
#---Figure 4
################################################################################
rm(list=ls())
# 加载需要的R包
library(Seurat)
library(monocle)
library(dplyr)
#Extract data, phenotype data, and feature data from the SeuratObject
load("~/Downloads/scRNA/snRNA.RData")
table(snRNA@meta.data$celltype)
AB <- snRNA[,snRNA@meta.data$celltype %in% c("NSC/Astrocyte","ImN","DG_ExN","InN")]
save(AB,file="/monocle/RGL_ImN_GC_InN.RData")
AB <- subset(AB, downsample = 3000)
data <- as(as.matrix(AB@assays[["RNA"]]@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = AB@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
cds  <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
HSMM <- cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
print(head(pData(HSMM)))
#select genes
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 10)
deg <- subset(diff_test_res, qval < 0.01)
deg <- deg[order(deg$qval, decreasing = F),]
head(deg)
write.table(deg, file="train.monocle.DEG.xls", col.names = T, row.names = F, sep = "\t", quote = F)
ordering_genes <- row.names (deg) 
HSMM <- setOrderingFilter(HSMM, ordering_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(HSMM)
dev.off()
#reduceDimension
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
#order
HSMM <- orderCells(HSMM)
HSMM@auxOrderingData[[HSMM@dim_reduce_type]]$branch_points
HSMM=orderCells(HSMM,root_state = 3) 

#plot
################################################################################
#---Figure 4A
################################################################################
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(HSMM,color_by = "celltype",cell_size = 0.1)
ggsave("celltype_all.pdf",device = "pdf",width = 10,height = 12,units = c("cm"))
plot_cell_trajectory(HSMM,color_by = "celltype",cell_size = 0.1)+facet_wrap(~Age,nrow=1)
ggsave("celltype_Age.pdf",device = "pdf",width = 22,height = 10,units = c("cm"))

################################################################################
#---Figure 4B
################################################################################
library(ggpubr)
df <- pData(HSMM)
ggdensity(df, x="Pseudotime",y="..density..", combine = TRUE,
          xlab = "Pseudotime", add = "median", rug = TRUE, color = "celltype", 
          fill = "celltype", palette = "jco",ncol=3)
ggsave("densirt.pdf",width = 10,height = 10,units = "cm")

################################################################################
#---Figure 4C
################################################################################
plot_cell_trajectory(HSMM, color_by = "orig.ident",cell_size = 0.2,use_color_gradient = TRUE,markers = c("SOX2","GFAP","PROX1","SNAP25","SLC1A3","PAX6","HOPX","PDGFRA","SNAP25","OLIG1","VCAN"))
ggsave("Marker_state.pdf",device = "pdf",width = 20,height = 20,units = c("cm"))

################################################################################
#---Figure 4D
################################################################################
count_table <- table(HSMM$age, HSMM$State)
count_table
#bar plot using Prism 9
plot_cell_trajectory(HSMM, color_by = "State",cell_size = 0.2)
ggsave("State.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))

################################################################################
#---Figure 4E
################################################################################
HSMM <- readRDS("HSMM.rds")
BEAM_res <- BEAM(HSMM,branch_point = 1,cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")
library('RColorBrewer')
tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4,
                                 cores = 1,
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T
)
pdf("branched_heatmap.pdf",width = 6,height = 6)
tmp1$ph_res
dev.off()

#GO enrichment analysis
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
write.csv(gene_group,"GP_gene_group.csv") 
write.csv(allcluster_go,"allcluster_go.csv") 
write.csv(go_res,"go_res.csv") 
write.csv(diff_test_res,"diff_test_res.csv")
saveRDS(HSMM, file = "HSMM.rds")

#cluster gene score
load("~/Downloads/snRNA/RGL_ImN_GC_InN.RData")
cluster_gene <- read_xlsx("cluster_gene.xlsx")
gene <- as.list(cluster_gene)
AB <- AddModuleScore(AB, features = gene, ctrl = 100, name = "cluster_gene")
colnames(AB@meta.data)[11]<-"cluster_gene_Score" 
#VlnPlot
P1<-VlnPlot(AB, features = 'cluster_Score',pt.size=0.01)
ggsave(filename = "cluster.pdf", plot = P1, device = 'pdf', width = 12, height = 10, units = 'cm')
rm('P1') 

