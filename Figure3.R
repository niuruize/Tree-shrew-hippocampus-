################################################################################
# Figure 3B
################################################################################
# Load R libraries
library(Seurat)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)
library(hopach)
  
# Load Aging Seurat object
load("~/Downloads/scRNA/snRNA.RData")
AB = snRNA
table(AB@meta.data[["celltype"]])
AB$grouping=str_replace(AB$Age,"_.*$","")
celltypes <- unique(AB@meta.data$celltype)
celltypes <- celltypes[which(!is.na(celltypes))]
celltypes <- setdiff(celltypes, c("IPC1","IPC2","Microglia","OPC","NFOL","MOL",
                                  "ImN","DG_ExN","nonDG_ExN","InN","CR",
                                  "End","Ependymal","Pericyte","Unk"))
getEuclideanDistance <- function(celltype, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp= AB[,AB@meta.data$celltype %in% c("NSC/Astrocyte")]
  expr <- tmp@assays[["RNA"]]@counts
  
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data$grouping)[c("Infancy", "Adult","Old")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Old")], nsample)
  adult_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Adult")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Infancy")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, adult_r,old_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, young, adult, old){
    tmp <- data.matrix(sqrt(matr[genes, young]))
    mean <- rowMeans(sqrt(matr[genes, young]))
    d_young <- distancevector(t(tmp), mean , d="euclid")
    names(d_young) <- young
    tmp <- data.matrix(sqrt(matr[genes, adult]))
    mean <- rowMeans(sqrt(matr[genes, adult]))
    d_adult <- distancevector(t(tmp), mean , d="euclid")
    names(d_adult) <- adult
    tmp <- data.matrix(sqrt(matr[genes, old]))
    mean <- rowMeans(sqrt(matr[genes, old]))
    d_old <- distancevector(t(tmp), mean , d="euclid")
    names(d_old) <- old
    
    list(young = d_young, adult=d_adult, old = d_old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, old = old_r, adult = adult_r, young = young_r)
  ds
}

# Run for all celltypes ####
res <- lapply(celltypes, function(x) getEuclideanDistance(x, lowcv = F))
names(res) <- celltypes
NSC_Astrocyte <- res
res<-c(NSC_Astrocyte,IPC1,IPC2,Microglia,OPC,NFOL,MOL,ImN,DG_ExN,nonDG_ExN,InN,CR,End,Ependymal,Pericyte,Unk)

res_original <- res
# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"
diffs_2 <- unlist(lapply(res_original, function(x) log2(mean(x[[3]]) / mean(x[[2]]))))
pvals_2 <- unlist(lapply(res_original, function(x) wilcox.test(x[[2]], x[[3]])$p.value))
adj_pvals_2 <- p.adjust(pvals_2, method = "BH")
save(res,file="res.RData")
pdf("Transcriptional_noise.pdf",width = 6,height = 5)
ord <- rev(order(diffs))
par(mar = c(15,5,2,5))
boxplot(do.call(c, res[ord]), las = 2, outline = F, col = c("blue", "red","green"), ylab = "Transcriptional noise", xaxt = 'n')
axis(1, at = seq(from = 1.5, to = 39.5, by = 3), names(res)[ord], las = 2)
dev.off()

################################################################################
# Figure 3D,E: Aging vs Adult in each cell type
################################################################################

library(MAST);library(Seurat);library(dplyr)
#SingleCellAssay object
Idents(snRNA)="celltype"
table(Idents(snRNA)) 
#DEGs were calculated separately for each cell type
snRNA1 = snRNA[,snRNA$celltype %in% c("ExN")]
fData = data.frame(symbolid=rownames(snRNA1),primerid=rownames(snRNA1))
rownames(fData)=fData$symbolid
cData = snRNA1@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(snRNA1@assays$RNA@data), cData = cData,fData = fData)
rm(snRNA1)
###------------------------Aging vs Adult--------------------###
gc()
dim(sca)
table(colData(sca)$Group)
cond<-factor(colData(sca)$Group)
cond<-relevel(cond,"Adult")
colData(sca)$condition<-cond

## (1) cngeneson
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + Sex, sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionAging')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

#
df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionAging') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionAging') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"ExN_MAST_Aging_Adult_DEGs.csv")

#the heatmap(D) and (E) were drawn using TBtools

################################################################################
# Figure 3H
################################################################################
#Aging_Up_score
load("~/Downloads/snRNA/snRNA.RData")
Aging_Up_gene <- read_xlsx("Aging_Up.xlsx")
gene <- as.list(Aging_Up_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Aging")
colnames(AB@meta.data)[11]<-"Human_aging_up_Score" 
#RidgePlot
P1<-RidgePlot(AB, features = 'Human_aging_up_Score', ncol = 1) 
ggsave(filename = "Human_aging_up_Ridge.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
#boxplot
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggboxplot(b, x = "orig.ident", y = "Human_aging_up_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "orig.ident", palette = "npg")+stat_compare_means(comparisons = my_comparisons)
ggsave("Human_aging_up_boxplot.pdf",width = 15,height = 15,units = "cm")

################################################################################
# Figure 3I
################################################################################
#Aging_Down_score
load("~/Downloads/snRNA/snRNA.RData")
Aging_Down_gene <- read_xlsx("Aging_Down.xlsx")
gene <- as.list(Aging_Down_gene)
AB <- AddModuleScore(snRNA, features = gene, ctrl = 100, name = "Aging_Down")
colnames(AB@meta.data)[11]<-"Human_aging_down_Score" 
#RidgePlot
P1<-RidgePlot(AB, features = 'Human_aging_down_Score', ncol = 1) 
ggsave(filename = "Human_aging_down_Ridge.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1')  
#boxplot
my_comparisons <- list( c("Infancy", "Adult"), c("Infancy", "Old"), c("Adult", "Old") )
ggboxplot(b, x = "orig.ident", y = "Human_aging_down_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "orig.ident", palette = "npg")+stat_compare_means(comparisons = my_comparisons)
ggsave("Human_aging_down_boxplot.pdf",width = 15,height = 15,units = "cm")

