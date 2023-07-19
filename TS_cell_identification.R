#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

#loading data
load("~/AB.RData")

################################################################################
#FindAllMarkers
################################################################################

markers <- FindAllMarkers(AB, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
write.table(markers_df,file="markers500.txt",quote=F,sep="\t",row.names=F,col.names=T)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file="markers10.txt",quote=F,sep="\t",row.names=F,col.names=T)
Heatmap <- DoHeatmap(AB,features = top10$gene) 
ggsave(filename = "Heatmap.pdf", plot = Heatmap, device = 'pdf', width = 100, height = 80, units = 'cm')
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
Heatmap <- DoHeatmap(AB,features = top5$gene) 
ggsave(filename = "Heatmap5.pdf", plot = Heatmap, device = 'pdf', width = 100, height = 50, units = 'cm')

################################################################################
#VlnPlot for marker gene
################################################################################

DefaultAssay(AB) <- "RNA"
Astro_marker<-VlnPlot(AB, features = c("AQP4","GFAP","GPR98","MASS1","SOX9","SLC1A2",	"SLC1A3",	"GLUL",	"VIM","GPC5","RYR3"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "Astro_marker.pdf", plot = Astro_marker, device = 'pdf', width = 15, height = 10, units = 'cm')
rm('Astro_marker')                    

Microglia_marker<-VlnPlot(AB, features = c("C1QA","CSF1R","P2RY12","C1QB","CD74","LRMDA","APBB1IP",	"PTPRC"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "Microglia_marker.pdf", plot = Microglia_marker, device = 'pdf', width = 15, height = 8, units = 'cm')
rm('Microglia_marker')  

Ex_marker<-VlnPlot(AB, features = c("RGS4","PROX1","SATB2","CAMK2A","SNAP25","LDB2"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "Ex_marker.pdf", plot = Ex_marker, device = 'pdf', width = 15, height = 6, units = 'cm')
rm('Ex_marker')  

Neuron_marker<-VlnPlot(AB, features = c("STMN2","SNAP25","RBFOX3"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "Neuron_marker.pdf", plot = Neuron_marker, device = 'pdf', width = 15, height = 6, units = 'cm')
rm('Neuron_marker') 

InN_marker<-VlnPlot(AB, features = c("GAD1","GAD2","SNAP25","SLC6A1","SLC32A1",	"LHFPL3",	 "PCDH15",	"SST","CALB2","LHX6","NR2F2"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "InN_marker.pdf", plot = InN_marker, device = 'pdf', width = 15, height = 10, units = 'cm')
rm('InN_marker')  

End_marker<-VlnPlot(AB, features = c("FLT1","CLDN5","RGS5","EBF1"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "End_marker.pdf", plot = End_marker, device = 'pdf', width = 15, height = 6, units = 'cm')
rm('End_marker') 

MOL_marker<-VlnPlot(AB, features = c("MOG","PLP1","ST18"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "MOL_marker.pdf", plot = MOL_marker, device = 'pdf', width = 20, height = 6, units = 'cm')
rm('MOL_marker')  

OPC_marker<-VlnPlot(AB, features = c("VCAN","PDGFRA",	"PCDH15",	"OLIG1",	"LHFPL3"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "OPC_marker.pdf", plot = OPC_marker, device = 'pdf', width = 15, height = 8, units = 'cm')
rm('OPC_marker')  

NFOL_marker<-VlnPlot(AB, features = c("FYN","BCAS1"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "NFOL_marker.pdf", plot = NFOL_marker, device = 'pdf', width = 15, height = 4, units = 'cm')
rm('NFOL_marker')  

RGL_marker<-VlnPlot(AB, features = c("SLC1A3","SOX2",	"GFAP","SOX6","SOX5","PAX6","HOPX","NEUROD1","SOX4","SOX11"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "RGL_marker.pdf", plot = RGL_marker, device = 'pdf', width = 15, height = 6, units = 'cm')
rm('RGL_marker')  

IPC_marker<-VlnPlot(AB, features = c("SLC1A3","SOX2",	"GFAP","SOX6","ASCL1","RFC4"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "IPC_marker.pdf", plot = IPC_marker, device = 'pdf', width = 15, height = 8, units = 'cm')
rm('IPC_marker')  

GC_marker<-VlnPlot(AB, features = c("PROX1","DCX","NIFK","NPY1R"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "GC_marker.pdf", plot = GC_marker, device = 'pdf', width = 15, height = 6, units = 'cm')
rm('GC_marker') 

RC_marker<-VlnPlot(AB, features = c("SNAP25","SLC6A1","CCK","LAMP5","SV2C","CNR1","SST","CALB2","PVALB","DCX","PROX1"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "RC_marker.pdf", plot = RC_marker, device = 'pdf', width = 20, height = 8, units = 'cm')
rm('RC_marker') 

all_marker<-VlnPlot(AB, features = c("SOX6","SOX2","SLC1A3","GFAP","PAX6","HOPX","AQP4","PTPRC","C1QA","C1QB","CSF1R","VCAN","PDGFRA",
                                     "MOG","PLP1","ST18","SNAP25","PROX1","DCX","SATB2","CAMK2A","SLC6A1","CCK","LAMP5","SV2C","CNR1",
                                     "SST","CALB2","PVALB","RELN","EBF1","RGS5"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "all_marker.pdf", plot = all_marker, device = 'pdf', width = 20, height = 30, units = 'cm')
rm('all_marker')

################################################################################
#Specifies the cell type for the cluster
################################################################################

#Doublets and cells with poor quality were removed
AB<- AB[,AB@meta.data$seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9","10",
                                             "11","12","13","14","15","16","17","18","19","20",
                                             "21","23","24","25","26","27","28","30",
                                             "31","32","33","34")]
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10",
                         "11","12","13","14","15","16","17","18","19","20",
                         "21","23","24","25","26","27","28","30",
                         "31","32","33","34")
new.cluster.ids <- c("MOL","NSC/Astrocyte","nonDG_ExN","DG_ExN","OPC","nonDG_ExN","Unk","nonDG_ExN","InN","InN","nonDG_ExN", #0-10
                     "DG_ExN","Microglia","DG_ExN","DG_ExN","nonDG_ExN","InN","InN","End","ImN","InN", #11-20
                     "nonDG_ExN","nonDG_ExN","InN","nonDG_ExN","NFOL","Ependymal","IPC1","nonDG_ExN","Pericyte","CR","ImN","IPC2") 
AB@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AB@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
Idents(AB) <- AB@meta.data$celltype
AB@meta.data[["celltype"]]<-factor(Ab@meta.data[["celltype"]], levels=c("NSC/Astrocyte","IPC1","IPC2","Microglia","OPC","NFOL","MOL",
                                          "ImN","DG_ExN","nonDG_ExN","InN","CR","End","Ependymal","Pericyte","Unk"))
AB@meta.data[["Age"]]<-factor(AB@meta.data[["Age"]], levels=c("Infancy","Adult","Old"))

