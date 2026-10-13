library(devtools)
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(hdf5r)
library(remote)
library(BiocManager)

# Load custom functions
source("src", "utils.R")

# Set working directory
data_dir <- file.path("data")
out_dir <- file.path("outs")

# Load the dataset
GBS1 <- Read10X_h5(file.path(data_dir, "1", "output_filtered.h5"), use.names = T)
GBS1 <- CreateSeuratObject(counts = GBS1, project = "pbmc", min.cells = 3, min.features = 200)
GBS1$condition <- "GBS"
GBS1$sample <- "1"
GBS1 <- NormalizeData(GBS1, verbose = FALSE)
GBS1 <- FindVariableFeatures(GBS1, selection.method = "vst", nfeatures = 2000)
GBS1[["percent.mt"]] <- PercentageFeatureSet(GBS1, pattern = "^MT-")
VlnPlot(GBS1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS1_QC_before filter", Sys.Date())), w = 10, h = 5)
GBS1 <- subset(GBS1, subset = 
                 nFeature_RNA > 800 & 
                 nFeature_RNA < 9000 & 
                 nCount_RNA < 30000 & 
                 percent.mt < 5)
VlnPlot(GBS1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS1_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(GBS1)[grepl("^MT-", rownames(GBS1))]
GBS1[["RNA"]] <- subset(GBS1[["RNA"]], features=setdiff(rownames(GBS1), mt))

C1 <- Read10X_h5(file.path(data_dir, "2", "output_filtered.h5"), use.names = T)
C1 <- CreateSeuratObject(counts = C1, project = "pbmc", min.cells = 3, min.features = 200)
C1$condition <- "Control" 
C1$sample <- "2"
C1 <- NormalizeData(C1, verbose = FALSE)
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C1_QC_before filter", Sys.Date())), w = 10, h = 5)
C1 <- subset(C1, subset = 
               nFeature_RNA > 200 & 
               nFeature_RNA < 6000 & 
               nCount_RNA < 15000 & 
               percent.mt < 5)
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C1_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(C1)[grepl("^MT-", rownames(C1))]
C1[["RNA"]] <- subset(C1[["RNA"]], features=setdiff(rownames(C1), mt))

GBS2 <- Read10X_h5(file.path(data_dir, "3", "output_filtered.h5"), use.names = T)
GBS2 <- CreateSeuratObject(counts = GBS2, project = "pbmc", min.cells = 3, min.features = 200)
GBS2$condition <- "GBS"
GBS2$sample <- "3"
GBS2 <- NormalizeData(GBS2, verbose = FALSE)
GBS2 <- FindVariableFeatures(GBS2, selection.method = "vst", nfeatures = 2000)
GBS2[["percent.mt"]] <- PercentageFeatureSet(GBS2, pattern = "^MT-")
VlnPlot(GBS2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS2_QC_before filter", Sys.Date())), w = 10, h = 5)
GBS2 <- subset(GBS2, subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 6000 & 
                 nCount_RNA < 20000 & 
                 percent.mt < 5)
VlnPlot(GBS2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS2_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(GBS2)[grepl("^MT-", rownames(GBS2))]
GBS2[["RNA"]] <- subset(GBS2[["RNA"]], features=setdiff(rownames(GBS2), mt))

C2 <- Read10X_h5(file.path(data_dir, "4", "output_filtered.h5"), use.names = T)
C2 <- CreateSeuratObject(counts = C2, project = "pbmc", min.cells = 3, min.features = 200)
C2$condition <- "Control"
C2$sample <- "4"
C2 <- NormalizeData(C2, verbose = FALSE)
C2[["percent.mt"]] <- PercentageFeatureSet(C2, pattern = "^MT-")
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C2_QC_before filter", Sys.Date())), w = 10, h = 5)
C2 <- subset(C2, subset = 
               nFeature_RNA > 200 & 
               nFeature_RNA < 7000 & 
               nCount_RNA < 20000 & 
               percent.mt < 7)
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C2_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(C2)[grepl("^MT-", rownames(C2))]
C2[["RNA"]] <- subset(C2[["RNA"]], features=setdiff(rownames(C2), mt))

C3 <- Read10X_h5(file.path(data_dir, "5", "output_filtered.h5"), use.names = T)
C3 <- CreateSeuratObject(counts = C3, project = "pbmc", min.cells = 3, min.features = 200)
C3$condition <- "Control"
C3$sample <- "5"
C3 <- NormalizeData(C3, verbose = FALSE)
C3[["percent.mt"]] <- PercentageFeatureSet(C3, pattern = "^MT-")
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C3_QC_before filter", Sys.Date())), w = 10, h = 5)
C3 <- subset(C3, subset = 
               nFeature_RNA > 200 & 
               nFeature_RNA < 6000 & 
               nCount_RNA < 20000 & 
               percent.mt < 5)
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("C3_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(C3)[grepl("^MT-", rownames(C3))]
C3[["RNA"]] <- subset(C3[["RNA"]], features=setdiff(rownames(C3), mt))

GBS3 <- Read10X_h5(file.path(data_dir, "7", "output_filtered.h5"), use.names = T)
GBS3 <- CreateSeuratObject(counts = GBS3, project = "pbmc", min.cells = 3, min.features = 200)
GBS3$condition <- "GBS" # Ish set GBS2$Condition to GBS here
GBS3$sample <- "7"
GBS3 <- NormalizeData(GBS3, verbose = FALSE)
GBS3[["percent.mt"]] <- PercentageFeatureSet(GBS3, pattern = "^MT-")
VlnPlot(GBS3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS3_QC_before filter", Sys.Date())), w = 10, h = 5)
GBS3 <- subset(GBS3, subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 6000 & 
                 nCount_RNA < 20000 & 
                 percent.mt < 7)
VlnPlot(GBS3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS3_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(GBS3)[grepl("^MT-", rownames(GBS3))]
GBS3[["RNA"]] <- subset(GBS3[["RNA"]], features=setdiff(rownames(GBS3), mt))

GBS4 <- Read10X_h5(file.path(data_dir, "8", "output_filtered.h5"), use.names = T)
GBS4 <- CreateSeuratObject(counts = GBS4, project = "pbmc", min.cells = 3, min.features = 200)
GBS4$condition <- "GBS"
GBS4$sample <- "8"
GBS4 <- NormalizeData(GBS4, verbose = FALSE)
GBS4[["percent.mt"]] <- PercentageFeatureSet(GBS4, pattern = "^MT-")
VlnPlot(GBS4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS4_QC_before filter", Sys.Date())), w = 10, h = 5)
GBS4 <- subset(GBS4, subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 6000 & 
                 nCount_RNA < 20000 & 
                 percent.mt < 4)
VlnPlot(GBS4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(path = file.path(out_dir, paste0("GBS4_QC_after filter", Sys.Date())), w = 10, h = 5)
mt <- rownames(GBS4)[grepl("^MT-", rownames(GBS4))]
GBS4[["RNA"]] <- subset(GBS4[["RNA"]], features=setdiff(rownames(GBS4), mt))

#perform integration
anchors <- FindIntegrationAnchors(object.list = list(GBS1,GBS2,GBS3,GBS4,C1,C2,C3), dims = 1:10)
combined <- IntegrateData(anchorset = anchors, dims = 1:10)

saveRDS(combined, file = file.path("outs", "GBS.rds"))

#Perform an integrated analysis
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
ElbowPlot(combined, ndims = 50)
saveThePlot(path = file.path(out_dir, paste0("Elbow", Sys.Date())), w = 7, h = 6)
DimHeatmap(combined, dims = 1:20, cells = 500, balanced = TRUE)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

# UMAP and Clustering

combined <- RunUMAP(combined, reduction = "pca", dims = 1:25)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:25)
combined <- FindClusters(combined, resolution = 0.6)

saveRDS(combined, file = file.path("outs", "GBS.rds"))