# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(rjson)
library(cowplot)

setwd("~/2025_KOGO_workshop/spatial/")


# Chapter 6. Xenium in situ --------


## Load the Xenium data --------
# xenium.obj <- LoadXenium('data/xenium/', fov = 'fov')
# saveRDS(xenium.obj, file = 'outs/xenium/xenium_obj.rds')
xenium.obj <- readRDS('outs/xenium/xenium_obj.rds')

## Preprocessing --------
# Remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# Add metadata
xenium.obj@meta.data$cells <- 'cells'


## Visualize the position --------
ImageDimPlot(xenium.obj, fov = "fov", 
             molecules = c("REG1A", "SPP1", "STAB1", "TGFBI"), 
             group.by = 'cells', nmols = 20000) 


## Visualize the expression level --------
ImageFeaturePlot(xenium.obj, 
                 features = c("CEACAM6", "PIGR"), 
                 max.cutoff = c(20, 20), size = 1.0, cols = c("white", "red"))



## Zoom in on a chosen area --------
# Increase your RAM usage (8GB)
options(future.globals.maxSize = 8000 * 1024^2)

# Define cropped area
ImageDimPlot(xenium.obj, fov = "fov", 
             molecules = c("REG1A", "SPP1", "STAB1", "TGFBI"), 
             axes = T, group.by = 'cells', nmols = 20000) 

cropped.coords <- Crop(xenium.obj[["fov"]], 
                       x = c(3000, 4000), 
                       y = c(2800, 3200), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords

# Visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"

ImageDimPlot(xenium.obj, fov = "zoom", 
             axes = TRUE, border.color = "white", border.size = 0.1, 
             cols = "polychrome", coord.fixed = FALSE, 
             molecules = c("REG1A", "SPP1", "STAB1", "TGFBI"), 
             nmols = 10000, group.by = 'cells')


## Preprocessing --------
# DO NOT RUN #
# xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
# xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
# xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
# xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)
# xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
# saveRDS(xenium.obj, file = 'outs/xenium/xenium_preprocessed.rds')

# Load preprocessed data
xenium.obj <- readRDS('outs/xenium/xenium_preprocessed.rds')


## Visualization --------
DimPlot(xenium.obj, reduction="umap")

FeaturePlot(xenium.obj, features = c("REG1A", "SPP1", "STAB1", "TGFBI"))

ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)

