# 19th KOGO single-cell RNA-seq Seurat practice #

# 0. Working directory Setting ----
getwd()
setwd('~/2025_KOGO_workshop/sc_rnaseq')

# 1. Load required Packages ----
install.packages('cowplot')
library(cowplot)
BiocManager::install('dittoSeq')
library(dittoSeq)
library(Seurat)
library(Azimuth)
library(BPCells)
library(ggplot2)
library(stringr)

# 2. Load data ----
## 1) count matrix ----
## Read Aggr'd dataset as a Seurat V5 object
FlexOutPath <- './data/' # Path to cellranger aggr output folder
ColonFlex.data <- open_matrix_10x_hdf5(path = paste0(FlexOutPath,'HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5')) 
write_matrix_dir(mat = ColonFlex.data,dir = './data/Outputs/FlexSeurat/') 
Flex.mat <- open_matrix_dir(dir = './data/Outputs/FlexSeurat/')

## Convert Ensemble name to gene symbol
Flex.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = Flex.mat, species = "human") 

## 2) metadata ----
# Read aggregation.csv file to be used as MetaData (Patient, etc)
MetaData <- read.csv(paste0(FlexOutPath,'HumanColonCancer_Flex_Multiplex_aggregation.csv'))
MetaData$Patient <- sapply(strsplit(MetaData$sample_id,"_"),function(X){return(X[6])})
MetaData$Patient_id <- str_sub(sapply(strsplit(MetaData$sample_id,"_"),function(X){return(X[6])}),1,2)
MetaData$Condition <- str_sub(sapply(strsplit(MetaData$sample_id,"_"),function(X){return(X[6])}),-3,-1)
Index <- as.numeric(sapply(strsplit(colnames(Flex.mat),"-"),function(X){return(X[2])}))
MetaData <- MetaData[Index,3:5]
MetaData$Barcode <- colnames(Flex.mat)
rownames(MetaData) <- MetaData$Barcode
MetaData %>% head

# 3. Work witha a Seurat object ----
## 1) Create Seurat Object ----
ColonCancer_Flex <- CreateSeuratObject(Flex.mat,meta.data = MetaData, project = "CRC")
ColonCancer_Flex
# An object of class Seurat 
# 18042 features across 279609 samples within 1 assay 
# Active assay: RNA (18042 features, 0 variable features)
# 1 layer present: counts

ColonCancer_Flex@assays$RNA$counts[1:10,1:5] %>% as.matrix()

## 2) Visualize the number of cell counts per sample ----
Idents(ColonCancer_Flex) <- ColonCancer_Flex$Patient
ColonCancer_Flex@meta.data %>% 
  ggplot(aes(x=Patient, fill=Patient)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = after_stat(count)),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")

saveRDS(ColonCancer_Flex, file = './results/ColonCancer_Flex.rds')

# 4. Random subsampling ----
## 1) Random subsampling on each sample ----
ColonCancer_Flex

set.seed(123) # to ensure that the same cells are subsampled

subset_cells <- ColonCancer_Flex@meta.data %>%
  group_by(ident = Idents(ColonCancer_Flex)) %>%
  sample_frac(size = 0.05) %>%
  pull(Barcode)

crc <- subset(ColonCancer_Flex, cells = subset_cells)

crc

table(crc$Patient)
# P1CRC P2CRC P2NAT P3CRC P3NAT P4CRC P5CRC P5NAT 
# 1025  2097  1640  2533  2062  1550  1231  1844 

## 2) Visualize the number of cell counts per sample ----
Idents(crc) <- crc$orig.ident
crc@meta.data %>% 
  ggplot(aes(x=Patient, fill=Patient)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = after_stat(count)),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")

saveRDS(crc, file = './results/crc.rds')

# 5. Standard pre-processing ----
crc <- readRDS('./results/crc.rds')

## 1) mitochondria counts ----
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
crc[["percent.mt"]] <- PercentageFeatureSet(crc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(crc@meta.data, 5)
# orig.ident nCount_RNA nFeature_RNA Patient Patient_id Condition                    Barcode percent.mt
# AAACAAGCACATAGTGACTTTAGG-1        CRC      13746         5165   P2CRC         P2       CRC AAACAAGCACATAGTGACTTTAGG-1  2.8080896
# AAACGGGCACTAAGGCACTTTAGG-1        CRC       1683         1264   P2CRC         P2       CRC AAACGGGCACTAAGGCACTTTAGG-1  0.5347594
# AAACTGGGTCTACCAGACTTTAGG-1        CRC        421          259   P2CRC         P2       CRC AAACTGGGTCTACCAGACTTTAGG-1  1.1876485
# AAAGCATGTTGTGAGAACTTTAGG-1        CRC       1915         1310   P2CRC         P2       CRC AAAGCATGTTGTGAGAACTTTAGG-1  0.5221932
# AAAGGGATCATTGTACACTTTAGG-1        CRC       1497          907   P2CRC         P2       CRC AAAGGGATCATTGTACACTTTAGG-1  2.5384102

## 2) Visualize QC metrics as a violin plot ----
VlnPlot(crc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## 3) Visualize as a FeatureScatter plot ----
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(crc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(crc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## 4) Quality control ----
# UMI and Gene Threshold (by removing 5%)
UMI_TH<-quantile(crc$nCount_RNA,c(0.025,0.975))
Gene_TH<-quantile(crc$nFeature_RNA,c(0.025,0.975))

# Add variable with QC filter status
crc$QCFilter<-ifelse(crc$percent.mt < 25 & 
                       crc$nCount_RNA > UMI_TH[1] & crc$nCount_RNA < UMI_TH[2] & 
                       crc$nFeature_RNA > Gene_TH[1] & crc$nFeature_RNA < Gene_TH[2],"Keep","Remove")

# Remove cells that failed QC
crc <- subset(crc,cells=colnames(crc)[crc$QCFilter=="Keep"])

saveRDS(crc, file = './results/filtered_crc.rds')

# 6. Normalization and find variable features ----
## 1) Normalization ----
# crc <- readRDS("./results/filtered_crc.rds")
crc <- NormalizeData(crc, normalization.method = "LogNormalize", scale.factor = 10000)
crc@assays$RNA$data[1:10,1:5] %>% as.matrix()

## 2) Feature selection ----
## Identification of highly variable features (feature selection) ----
crc <- FindVariableFeatures(crc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(crc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(crc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# 7. Scaling the data ----
crc <- ScaleData(crc)
crc
crc@assays$RNA$scale.data[1:10,1:5] %>% as.matrix()

# 8. Cell cycle scoring ----
## 1) prepare cell cycle related genes ----
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## 2) do cell cycle scoring ----
crc <- CellCycleScoring(crc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(crc[[]])

## 3) Perform linear dimensional reduction ----
crc <- RunPCA(crc, features = VariableFeatures(object = crc))

## 4) Plot the PCA colored by cell cycle phase ----
DimPlot(crc,
        reduction = "pca",
        group.by= "Phase")

DimPlot(crc,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# 10. Determine the 'dimensionality' of the dataset ----
ElbowPlot(crc, ndims = 40)

# Determine percent of variation associated with each PC
pct <- crc[["pca"]]@stdev / sum(crc[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
pcs <- which(cumu > 90 & pct < 5)[1]
pcs # 42

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() + theme(aspect.ratio = 1)

# 11. Cluster the cells ----
## 1) Find neighbors ----
crc <- FindNeighbors(crc, dims = 1:42)

## 2) Find clusters ----
### 0.1, 0.2, 0.3, 0.4, 0.5, 0.6
crc <- FindClusters(crc, resolution =c(seq(0.1,0.6,0.1)))

## 3) Dimensional reduction using UMAP ----
crc <- RunUMAP(crc, dims = 1:42)

## 4) Visualize clustering results ----
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.1") + theme(aspect.ratio = 1)
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.2") + theme(aspect.ratio = 1)
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.3") + theme(aspect.ratio = 1)
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.4") + theme(aspect.ratio = 1)
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.5") + theme(aspect.ratio = 1)
DimPlot(crc, reduction = "umap", group.by = "RNA_snn_res.0.6") + theme(aspect.ratio = 1)

## 5) Visualize various effects ----
DimPlot(crc, reduction = "umap",
        group.by = c("RNA_snn_res.0.1","Phase","Patient_id","Condition"), ncol = 2) + theme(aspect.ratio = 1)

# 12. Finding differentially expressed features ----
## 1) Find markers for all clusters ----
Idents(crc) <- crc$RNA_snn_res.0.1
Mks <- FindAllMarkers(crc,min.diff.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
Mks %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1) %>% View()

## 2) make an expression heatmap (top 5 markers for each cluster)
Mks %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

data.frame(top5$cluster,top5$gene)
DoHeatmap(crc, features = top5$gene) + NoLegend()

## 3) Visualization using FeaturePlot and VlnPlot
Mks %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC > 1) %>%
  slice_head(n = 1) %>%
  ungroup() -> top1

top1

f1 <- FeaturePlot(crc, features = (top1[1,] %>% pull(gene))) + theme(aspect.ratio = 1)
v1 <- VlnPlot(crc, features = (top1[1,] %>% pull(gene))) + geom_boxplot(alpha=0.5) + theme(aspect.ratio = 1)
plot_grid(f1,v1, rel_widths = c(1,1.2))

f2 <- FeaturePlot(crc, features = (top1[2,] %>% pull(gene))) + theme(aspect.ratio = 1)
v2 <- VlnPlot(crc, features = (top1[2,] %>% pull(gene))) + geom_boxplot(alpha=0.5) + theme(aspect.ratio = 1)
plot_grid(f2,v2, rel_widths = c(1,1.2))

f3 <- FeaturePlot(crc, features = (top1[3,] %>% pull(gene))) + theme(aspect.ratio = 1)
v3 <- VlnPlot(crc, features = (top1[3,] %>% pull(gene))) + geom_boxplot(alpha=0.5) + theme(aspect.ratio = 1)
plot_grid(f3,v3, rel_widths = c(1,1.2))

f4 <- FeaturePlot(crc, features = (top1[4,] %>% pull(gene))) + theme(aspect.ratio = 1)
v4 <- VlnPlot(crc, features = (top1[4,] %>% pull(gene))) + geom_boxplot(alpha=0.5) + theme(aspect.ratio = 1)
plot_grid(f4,v4, rel_widths = c(1,1.2))

# 13. Rename Clusters (Manual Annotation) ----
crc <- RenameIdents(object = crc,
                    '0' = 'Tumor',
                    '1' = 'Smooth muscle 1',
                    '2' = 'Fibroblast',
                    '3' = 'T cells',
                    '4' = 'Intestinal epithelial',
                    '5' = 'Myeloid',
                    '6' = 'B cells 1',
                    '7' = 'B cells 2',
                    '8' = 'Endothelial 1',
                    '9' = 'Neuronal',
                    '10' = 'Myeloid',
                    '11' = 'Endothelial 2',
                    '12' = 'Smooth muscle 2')
DimPlot(crc, label = TRUE) + theme(aspect.ratio = 1)

# 14. Visualization of cluster proportion for each sample ----
crc$cluster_name <- Idents(crc)
dittoBarPlot(crc, "cluster_name", group.by = "Patient") + ggtitle("Cluster proportion by patient")

saveRDS(crc, file = './results/final_crc.rds')
