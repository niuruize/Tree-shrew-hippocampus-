################################################################################
# Supplementary Figure 2A
################################################################################
%%R
# too-many-cell-counts
snRNA <- subset(snRNA, downsample=1000)
count = as.data.frame(as.matrix(snRNA_F@assays$integrated@data))
write.csv(count,"count.csv")
label_celltype_species <- data.frame(item=colnames(count), label=snRNA$celltype_species)
write.csv(label_celltype_species,"labels.csv", row.names = F)

%%bash
docker run -it --rm -v /Users/niuruize/Downloads/TOO:/TOO \
      gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
    --matrix-path /TOO/input/count.csv \
    --labels-file /TOO/input/labels.csv \
    --smart-cutoff 1 \
    --min-size 1 \
    --draw-collection "PieRing" \
    --dendrogram-output "dendrogram.pdf" \
    --output /TOO/out \
    > clusters.csv

################################################################################
# Supplementary Figure 2B,C
################################################################################
%%R
# too-many-cell-counts
snRNA <- subset(snRNA, downsample=1000)
count = as.data.frame(as.matrix(snRNA_F@assays$integrated@data))
write.csv(count,"count.csv")
label_species <- data.frame(item=colnames(count), label=snRNA$species)
write.csv(label_species,"label_species.csv", row.names = F)

%%bash
docker run -it --rm -v /Users/niuruize/Downloads/TOO:/TOO \
      gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
    --matrix-path /TOO/input_species/count.csv \
    --labels-file /TOO/input_species/label_species.csv \
    --smart-cutoff 1 \
    --min-size 1 \
    --draw-collection "PieRing" \
    --dendrogram-output "dendrogram.pdf" \
    --output /TOO/out_species \
    > clusters.csv


