# load required packages
if(!require(pak)) install.packages("pak")
if(!require(Seurat)) remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
if(!require(BPCells)) devtools::install_github(c("bnprks/BPCells"))
library(BPCells)
if(!require(Azimuth)) remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(rstatix)
library(enrichR)
library(ggpubr)
if(!require(glmGamPoi)) BiocManager::install("glmGamPoi")
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
if(!require(ggvenn)) devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
if(!require(AUCell)) BiocManager::install("AUCell")
library(AUCell)
if(!require(pochi)) devtools::install_github("diegoalexespi/pochi")
library(pochi)
if(!require(WGCNA)) install.packages("WGCNA")
library(WGCNA)
library(scales)

#for Seurat.utils
if(!require(checkmate)) install.packages("checkmate")
if(!require(enrichplot)) BiocManager::install("enrichplot")

# install.packages("devtools"); # If you don't have it
if(!require(sparseMatrixStats)) BiocManager::install("sparseMatrixStats")
require("devtools")

# Install dependencies
if(!require(Stringendo)) devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
if(!require(CodeAndRoll2)) devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
if(!require(ReadWriter)) devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
if(!require(MarkdownHelpers)) devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
if(!require(Markdownreports)) devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
if(!require(ggExpress)) devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
if(!require(DatabaseLinke.R))devtools::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)
if(!require(Seurat.utils)) devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

# Import custom functions
source(file.path("src", "utils.R"))

# set paths
data_dir <- file.path("data") # data
dir <- file.path("outs", "Analysis_20260513", "Clusters17+18Remove") # outs

# load Seurat rds file
data <- LoadSeuratRds(file.path(dir, "GBS_Data_2024-11-15.rds"))

# load metadata with file paths, patient conditions, etc.
gbs_meta <- read.csv(file.path(data_dir, "gbs_metadata.csv"))

# QC
write.csv(table(data$sample), file=file.path(dir, "ncells_by_sample_original.csv"))

Idents(data) <- "sample"
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(umap, path = file.path(dir, paste0("QC_before filter_", Sys.Date())), w = 10, h = 5)
data <- subset(data, subset = 
                 nFeature_RNA > 200 & 
                 nFeature_RNA < 6000 & 
                 nCount_RNA < 25000 & 
                 percent.mt < 10)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveThePlot(umap, path = file.path(dir, paste0("QC_after filter_", Sys.Date())), w = 10, h = 5)
write.csv(table(data$sample), file=file.path(dir, "ncells_by_sample_filtered.csv"))

# visualise umaps
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", split.by = "sample", raster=F, label = T, ncol = 4)
saveThePlot(umap, path = file.path(dir, paste0("umap_by_sample", Sys.Date())), w = 20, h = 10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", split.by = "condition", raster=F, label = T, ncol = 4)
saveThePlot(umap, path = file.path(dir, paste0("umap_by_condition", Sys.Date())), w = 20, h = 10)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", raster=F, label = T)
saveThePlot(umap, path = file.path(dir, paste0("umap_unannotated", Sys.Date())), w = 9, h = 5)

# umaps with one cluster highlighted
for(i in 0:21) {
  cluster <- i
  cells.to.highlight <- names(Idents(data))[Idents(data)==cluster] # select all cells within the named cluster
  cluster.colour <- hue_pal()(length(levels(data@active.ident)))[cluster+1] # get cluster colour
  DimPlot(data, reduction = "umap", group.by = "seurat_clusters", 
          cells.highlight = cells.to.highlight,
          cols.highlight = cluster.colour,
          raster=F, label = T)
  saveThePlot(umap, path = file.path(dir, paste0("umap_unannotated_highlight_", 
                                                 cluster, "_", Sys.Date())), w = 8, h = 6)
}

# get # of cells per cluster per sample
n_cells <- table(Idents(data), data$sample) |> as.data.frame() |> 
  dplyr::rename("cluster" = Var1, "sample" = Var2, "n" = Freq) |>
  mutate(sample = as.numeric(sample)) |> 
  mutate(sample = ifelse(sample > 5, sample+1, sample)) |>
  left_join(gbs_meta, by = join_by(sample)) |>
  group_by(sample) |>
  mutate(prop = n/sum(n)*100)

n_cells$sample <- factor(n_cells$sample, levels = c("1", "3", "7", "8", "2", "4", "5"))

# plot n cells per cluster per sample
n_cells |> ggplot(aes(
  x = sample,
  y = prop,
  fill = cluster
)) + geom_col(position = "stack", color="black") +
  #scale_fill_brewer(palette = "Blues") +
  theme_classic() +
  scale_y_continuous(limits = c(0,101), expand = expansion(c(0, 0))) +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(colour = "white", fill = "white"))
saveThePlot(path = file.path(dir, paste0("abundance_by_sample")), w = 4, h = 6)

# plot n cells per cluster per CONDITION
n_cells <- table(Idents(data), data$condition) |> as.data.frame() |> 
  dplyr::rename("cluster" = Var1, "condition" = Var2, "n" = Freq) |>
  group_by(condition) |>
  mutate(prop = n/sum(n)*100)

n_cells |> ggplot(aes(
  x = condition,
  y = prop,
  fill = cluster
)) + geom_col(position = "stack", color="black") +
  theme_classic() +
  scale_y_continuous(limits = c(0,100)) +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(colour = "white", fill = "white"))
saveThePlot(path = file.path(dir, paste0("abundance_by_condition")), w = 2.5, h = 4)

# abundance plots
n_plot <- n_cells |> ggplot(aes(
  x = condition,
  y = prop
)) + stat_summary(aes(fill = condition), geom = "col", fun = "mean", position = "dodge", colour = "black") +
  scale_fill_brewer(palette = "Blues") +
  geom_point(colour = "black", size = 1,
             position = position_jitter(width = 0.2)) +
  scale_shape_manual(values=c(2,3)) +
  stat_summary(geom = "errorbar", width = 0.9, fun.data = mean_se, position = "dodge") + 
  facet_wrap(vars(cluster), scales = "free_y", axes = "all") + theme_classic() +
  stat_pvalue_manual(daa_t, label = "p", hide.ns = F, tip.length = 0.1) +
  scale_y_continuous(expand = expansion(c(0, 0.3))) +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(colour = "white", fill = "white"))
saveThePlot(path = file.path(dir, paste0("abundance_by_condition")), w = 7, h = 10)

##### annotate clusters #####
# identify clusters
# some panels to use for specific purposes
identity.general <- c("CD79A", "MS4A1", "CD3D", "CD4", "IL7R", "FOXP3", "CD8A", "GNLY", "NKG7", "NCAM1",
                      "CCL5", "SELL", "CCR7", "GIMAP5", "NRP1", "THBD", "CD1C", "CD14", "FCGR3A", 
                      "ITGB3", "ITGAX", "PECAM1", "MKI67", "TBX21", "GATA3", "RORC", 
                      "TNFRSF13B", "IGHD")

identity.tcell <- c("CD3G", "CD4", "CD8A", "CCR7", "SELL", "FAS", "HNRNPLL", "GIMAP5", 
                    "CD27", "CD28", "CD38", "CD69", 
                    "CXCR4", "CXCR5", "CCR4", "CCR6", "ARHGEF10", "PASK", "KLRB1", 
                    "TRAC", "TRBC2", "TRDC", "TRGC2",
                    "NLRP3", "CCL5", "CCL20", "TNF", "IL1B",
                    "TBX21", "IFNG", "GATA3", "RORC", "FOXP3", "IL7R", "IL2RA", "ID2", 
                    "TOX", "PRDM1", "BCL6", "PRF1", "NKG7", "GNLY", "GZMB", "GZMK", "GZMH",
                    "ICOS", "ITGAM", "HAVCR2", "PDCD1", "ENTPD1", "PECAM1", "NCAM1", 
                    "MKI67", "CDK1", "CDK4", "CCND2"
)

identity.cd4 <- c("CCR7", "SELL", "FAS", "HNRNPLL", "GIMAP5", 
                  "CD27", "CD28", "CD38", "CD44", "CD69", 
                  "CXCR4", "CXCR5", "CCR4", "CCR6", "ARHGEF10", "PASK", "KLRB1", 
                  "TRAC", "TRBC2", "TRDC", "TRGC2",
                  "NLRP3", "CCL5", "CCL20", "TNF", "IL1B",
                  "TBX21", "IFNG", "GATA3", "RORC", "FOXP3", "IL7R", "IL2RA", "ID2", 
                  "TOX", "PRDM1", "BCL6", "IL4",
                  "ICOS", "ITGAM", "HAVCR2", "PDCD1", "PECAM1", "NCAM1", 
                  "MKI67", "CDK1", "CDK4", "CCND2", "INPP4B"
)

identity.mono <- c(
  "CD14", "LYZ", "S100A8", "S100A9", "VCAN", 
  "FCGR3A", "TCF7L2", 
  "CXCR1", "CXCR2", "CLEC4E", "FCGR3B", "CSF3R",
  "GZMA", "ZAP70", "CD2", "CTSW", 
  "ITGB3"
)

identity.bcell <- c("IGKC", "IGLC1", #Ig kappa vs lambda light chain
                    "IGHD", # naive
                    "IGHM", # early memory
                    "TNFRSF13B", # memory
                    "CD27",
                    "SLAMF7", # plasma cell
                    "CD79A", "MS4A1"
)

identity.main <- c("CD3D", "CD4", "CD8A", "CCR7", "SELL", "FAS", "FOXP3",
                   "CCNL2", "PRF1", "GZMB", "GZMK", "GZMH", "GNLY", "MKI67", 
                   "NCAM1", "MS4A1", "CD79A", "IGHD", "TNFRSF13B", "SLAMF7", 
                   "CD14", "ITGAX", "IER2", "FCGR3A", "ITGB3", "THBD", "NRP1", "CD34"
)

#select which definition above to use
identity.markers <- identity.main

#dotplot
DotPlot(data, features = identity.markers, cluster.idents = T) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
saveThePlot(path = file.path(dir, paste0("feature_plot_dot_main_", Sys.Date())), w = 12, h = 7)

#heatmap
library("pheatmap")
hm.data <- AverageExpression(data)$RNA[identity.markers, ]
hm.data <- hm.data[,c("g0","g2","g11","g15","g7","g8","g1","g12",
                      "g18","g4","g5","g10","g20","g3",
                      "g6","g9","g13","g14","g16","g17","g19","g21")]

hm <- pheatmap(hm.data, cluster_rows=F, show_rownames=TRUE,
               cluster_cols=F, scale = "row",
               cellwidth = 9, cellheight = 9,
               filename = file.path(dir, paste0(
                 "feature_plot_pheatmap_main_", Sys.Date(), ".pdf")
               ))
hm
save_pheatmap_pdf(hm, file.path(dir, paste0("feature_plot_pheatmap_main_", Sys.Date(), ".pdf")),  
                  width = 6, height = 10)

# feature plot
feature_p <- FeaturePlot(data, features = c("OGFR", "IL6ST", "CLU", "GPR107"), raster = T, reduction = "umap")
feature_p
saveThePlot(feature_p, path = file.path(dir, paste0("feature_plot_umap_", Sys.Date())), w = 15, h = 25)

RidgePlot(data, features = identity.markers, ncol = 4)
saveThePlot(path = file.path(dir, paste0("ridgeline_plot_", Sys.Date())), w = 20, h = 40)


# find conserved markers across conditions in clusters and visualise
markers <- FindConservedMarkers(data, ident.1 = 1, grouping.var = "condition", verbose = TRUE)

conserved_p <- DotPlot(data, features = identity.markers, cols = c("blue", "red"), 
                       dot.scale = 8, split.by = "condition") +
  RotatedAxis()
conserved_p

# find positive markers for each cluster compared to all other remaining cells
markers <- FindAllMarkers(data, only.pos = TRUE)
write.csv(markers, file=file.path(dir, paste0("FindAllMarkers.csv")))
markers_filt <- markers |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)
write.csv(markers_filt, file=file.path(dir, paste0("FindAllMarkers_Filt.csv")))

cluster.to.enrich <- 0
cluster_enrich <- markers |> dplyr::filter(cluster==cluster.to.enrich & avg_log2FC > 1 & p_val_adj < 0.05) #|> head(n = 50)
cluster_enrich
write.csv(cluster_enrich, file=file.path(dir, paste0("Cluster", cluster.to.enrich, "_difffrom_all.csv")))


# find positive markers for a specific cluster compared to other specific clusters
find.marker.cluster <- c(1)
compare <- c(8,12)
find.marker.cluster_s <- paste(find.marker.cluster, collapse=",")
compare_s <- paste(compare, collapse=",")
markers_cluster <- FindMarkers(data, ident.1 = find.marker.cluster, 
                               ident.2 = compare) |> filter(avg_log2FC > 0.585 & p_val_adj < 0.05)
head(markers_cluster, n = 20)
write.csv(markers_cluster, file=file.path(dir, paste0("Cluster", 
                                                      find.marker.cluster_s, 
                                                      "_difffrom_",
                                                      compare_s,
                                                      ".csv")))
# volcano of DE genes for the cluster
# heatmap
markers.mono.filt <- markers.mono |> filter(avg_log2FC > 0.5 & p_val == 0) |>
  group_by(cluster) |> mutate(row = row_number(p_val_adj)) |> filter(row <=20)
DoHeatmap(subset(mono, downsample = 100), features = markers.mono.filt$gene, size = 3)
saveThePlot(path = file.path(dir, paste0("heatmap_mono_", Sys.Date())), w = 10, h = 10)

markers_cluster <- markers.mono |> filter(cluster==6)

# enrichr for cluster identity
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "KEGG_2016", "Panther_2016")
#enriched <- enrichr(as.vector(cluster_enrich$gene), dbs) # against all clusters from FindAllMarkers
enriched <- enrichr(as.vector(rownames(markers_cluster)), dbs) # for specific clusters from FindMarkers
write.csv(enriched[["GO_Molecular_Function_2023"]], file=file.path(dir, paste0("enrichment_GO_MF_", 
                                                                               find.marker.cluster_s, "_", compare_s, ".csv")))
write.csv(enriched[["GO_Biological_Process_2023"]], file=file.path(dir, paste0("enrichment_GO_BP_", 
                                                                               find.marker.cluster_s, "_", compare_s, ".csv")))
write.csv(enriched[["KEGG_2016"]], file=file.path(dir, paste0("enrichment_KEGG_", 
                                                              find.marker.cluster_s, "_", compare_s, ".csv")))
write.csv(enriched[["Panther_2016"]], file=file.path(dir, paste0("enrichment_Panther_", 
                                                                 find.marker.cluster_s, "_", compare_s, ".csv")))
plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("enrichment_GO_MF_", find.marker.cluster_s, "_", compare_s)), w = 8, h = 6)
plotEnrich(enriched[[2]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("enrichment_GO_BP_", find.marker.cluster_s, "_", compare_s)), w = 8, h = 6)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("enrichment_KEGG_", find.marker.cluster_s, "_", compare_s)), w = 8, h = 6)
plotEnrich(enriched[[4]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("enrichment_Panther_", find.marker.cluster_s, "_", compare_s)), w = 8, h = 6)

# save original Seurat object
saveRDS(data, file = file.path(dir, paste0("GBS_Data_", Sys.Date(), ".rds")))

# remove clusters 17 and 18
cells_to_remove <- WhichCells(data, idents = c(17, 18))
data <- subset(data, cells = setdiff(Cells(data), cells_to_remove))
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6",
                     "7", "8", "9", "10", "11", "12",
                     "13", "14", "15", "16", "17", "18", "19", "20",
                     "21") # overall


names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data[["seurat_clusters"]] <- Idents(data)
DimPlot(data, label = TRUE, pt.size = 0.5, raster=F, reduction = "umap") + NoLegend()
saveThePlot(path = file.path(dir, paste0("umap_annot_clean_", Sys.Date())), w = 9, h = 8)

# Broad cluster annotations
Idents(data) <- "seurat_clusters"
new.cluster.ids <- c("CD4 T", "CD8 T", "CD4 T", "Monocyte", "NK cell", "B cell", 
                     "Monocyte", "CD8 T", "CD8 T", "Monocyte", "B cell", "CD4 T",
                     "CD8 T", "Monocyte", "PMA", "CD4 T", "Platelet", "cDC", "NKT",
                     "pDC", "Plasma cell", "Progenitor")
new.cluster.ids <- c("0", "0", "1", "1", "4", "5", # CD4 T
                     "3", "1", "2", "2", "3", "2",
                     "3")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data[["broad.cluster"]] <- Idents(data)
DimPlot(data, label = TRUE, pt.size = 0.5, raster=F, reduction = "umap") + NoLegend()
saveThePlot(path = file.path(dir, paste0("umap_annot_broad_", Sys.Date())), w = 9, h = 8)

saveRDS(data, file = file.path(out_dir, "Analysis_20241007", "Clusters17+18Remove", 
                               paste0("GBS_Data_", Sys.Date(), ".rds")))


#### DIFFERENTIAL ANALYSIS #####
DefaultAssay(data) <- "RNA"

# differentially expressed genes
diff.all <- FindMarkers(data, 
                        ident.1 = "Control", 
                        ident.2 = "GBS",
                        test.use= "wilcox", 
                        min.pct=0.1, 
                        verbose = T,
                        only.pos=T)

# create pseudobulk for DEGs by condition (used for plotting the next step)
features.from.bulk_cd11b <- read.csv(file="11_sig_differential_expression_cd11b.csv") |> pull(symbol)
features.from.bulk_cd8 <- read.csv(file="11_sig_differential_expression_cd8.csv") |> pull(symbol)
features.from.bulk <- c(features.from.bulk_cd8, features.from.bulk_cd11b) |> unique()
features.from.bulk <- features.from.bulk[features.from.bulk %in% rownames(GetAssayData(data))] # limit to genes present in the sc data


# find differentially expressed genes per cluster and plot the biggest changes
top <- data.frame()
for(cluster in data$seurat_clusters |> levels()) {
  bulk <- subset(data, seurat_clusters == cluster)
  bulk.subset <- AggregateExpression(bulk, group.by = c("condition", "sample"), 
                                     slot = "counts", 
                                     assays = "RNA",  
                                     normalization.method = "LogNormalize",
                                     scale.factor = 10000,
                                     return.seurat = F)
  if(ncol(bulk.subset$RNA) < 3) next # skip iteration if  cluster is represented in less than 3 samples
  counts <- as.matrix(bulk.subset$RNA)
  coldata <- as.matrix(data.frame(condition = c("Control", "Control", "Control", "GBS", "GBS", "GBS", "GBS"), 
                                  row.names = c("Control_2", "Control_4", "Control_5", "GBS_1", 
                                                "GBS_3", "GBS_7", "GBS_8")))
  coldata <- subset(coldata, rownames(coldata) %in% colnames(counts)) # filters out samples where the cluster is not represented
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  
  # LFC SHRINKAGE
  resLFC <- lfcShrink(dds, coef=2, type="apeglm")
  resLFC
  
  plotMA(resLFC, ylim=c(-3,3), cex=.8)
  abline(h=c(-1,1), col="dodgerblue", lwd=2)
  
  diff <- as.data.frame(resLFC)
  diff
  
  diff["gene"] <- rownames(diff) # add gene names as their own column for easy referencing
  head(diff, n = 10)
  diff_write <- diff
  write.csv(diff_write, file=file.path(dir, paste0("diff_bulk_", cluster, "_result_allgenes.csv")))
  write.csv(diff_write |> filter(pvalue < 0.05), file=file.path(dir, paste0("diff_bulk_", cluster, "_result.csv")))
  #mark which genes are significantly different
  diff <- diff |>
    mutate(
      significance =
        case_when(
          padj < 0.05 & log2FoldChange > 0.5 ~ "increased",
          padj < 0.05 & log2FoldChange < -0.5 ~ "decreased",
          #TRUE ~ NA_character_,
          .default = "no change"
        )
    )
  
  # plot differentially expressed genes
  ggplot(diff, aes(log2FoldChange, -log10(pvalue), colour = significance)) + 
    geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
    
    ylab("-log10(adjusted p-value)") + 
    geom_text_repel(aes(label = ifelse(padj < 0.05, gene, ""))) + #, colour = "red", size = 3) +
    scale_colour_manual(
      values = c("decreased" = "blue", "increased" = "red", "no change" = "#cdcdcd"),
      na.value = "#cdcdcd"
    )
  saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_volcano")), w = 12, h = 8)
  
  #EnhancedVolcano
  diff.volc <- diff |> filter(!is.na(padj)) # filter out NAs caused by independent filtering of lowly expressed genes
  EnhancedVolcano(diff.volc,
                  lab = rownames(diff.volc),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  ylim = c(0, max(-log10(diff.volc[['pvalue']]), na.rm = TRUE) + 0.5),
                  xlim = c(-5, 8),
                  title = paste0('Cluster ', cluster),
                  subtitle = "",
                  drawConnectors = T,
                  pCutoff = 0.05,
                  pCutoffCol = 'padj',
                  FCcutoff = 0.5,
                  pointSize = 3.0,
                  labSize = 6.0)
  saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_volcano_enhanced")), w = 8, h = 8)
  
  
  violin_data <- as.data.frame(bulk.subset$RNA)
  violin_data["gene"] <- rownames(violin_data)
  violin_data <- violin_data |> pivot_longer(
    cols = starts_with("Control") | starts_with("GBS"),
    names_to = c("condition", "sample"),
    names_pattern = "(.*)_(.)",
    values_to = "count"
  )
  top_genes <- head(diff |> arrange(padj), n = 12)
  top_genes["gene"] <- rownames(top_genes)
  top_genes <- top_genes$gene
  violin_data |> filter(gene %in% top_genes) |> ggplot(aes(
    x = condition,
    y = count
  )) + geom_boxplot() + facet_wrap("gene", scales = "free")
  saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_boxplot")), w = 14, h = 8)
  
  # add to overall
  top <- rbind(top, head(diff |> arrange(padj), n = 10) |> mutate(cluster = cluster))
  
  # heatmap
  ntd <- normTransform(dds)
  #select <- order(rowMeans(counts(dds, normalized=TRUE)),
  #                decreasing=TRUE)[1:50]
  if(length(unique(diff$significance)) > 1) { # only generate heatmap if there are DE genes
    sel <- diff |> filter(significance != "no change") |> 
      arrange(desc(significance), padj) |> pull(gene)
    select <- match(sel, rownames(dds))
    if(length(select) > 1) { # if more than one DE gene
      df <- as.data.frame(colData(dds)["condition"]) |> arrange(condition)
      df2 <- data.frame(
        significance=diff$significance,
        row.names=diff$gene
      )
      library("pheatmap")
      hm <- pheatmap(assay(ntd)[select,], cluster_rows=F, show_rownames=TRUE,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                     cluster_cols=F, 
                     annotation_col=df, 
                     annotation_row=df2)
      hm
      save_pheatmap_pdf(hm, file.path(dir, paste0("diff_bulk_", cluster, "_Heatmap.pdf")),  
                        width = 6, height = 10)
    }
  }
  
  # enrichR for specifically this pseudobulked cluster
  markers_dea <- diff |> filter(log2FoldChange > 0.5 & padj < 0.05)
  if(dim(markers_dea)[1] > 0) { # check that df is not empty
    dbs <- listEnrichrDbs()
    dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "KEGG_2016", "Panther_2016")
    enriched <- enrichr(markers_dea$gene, dbs)
    enriched[["GO_Molecular_Function_2023"]]
    enriched[["GO_Biological_Process_2023"]]
    enriched[["KEGG_2016"]]
    enriched[["Panther_2016"]]
    if(dim(enriched[["GO_Molecular_Function_2023"]])[1] > 0) {
      plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
      saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_GO_MF_")), w = 8, h = 6)
      write.csv(enriched[[1]], file = file.path(dir, paste0("diff_bulk_", cluster, "_GO_MF.csv")))
    }
    if(dim(enriched[["GO_Biological_Process_2023"]])[1] > 0) {
      plotEnrich(enriched[[2]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
      saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_GO_BP_")), w = 8, h = 6)
      write.csv(enriched[[2]], file = file.path(dir, paste0("diff_bulk_", cluster, "_GO_BP.csv")))
    }
    if(dim(enriched[["KEGG_2016"]])[1] > 0) {
      plotEnrich(enriched[[3]], showTerms = 50, numChar = 60, y = "Ratio", orderBy = "P.value")
      saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_KEGG_")), w = 8, h = 6)
      write.csv(enriched[[3]], file = file.path(dir, paste0("diff_bulk_", cluster, "_KEGG.csv")))
    }
    if(dim(enriched[["Panther_2016"]])[1] > 0) {
      plotEnrich(enriched[[4]], showTerms = 5, numChar = 60, y = "Ratio", orderBy = "P.value")
      saveThePlot(path = file.path(dir, paste0("diff_bulk_", cluster, "_Panther_")), w = 8, h = 6)
      write.csv(enriched[[4]], file = file.path(dir, paste0("diff_bulk_", cluster, "_Panther.csv")))
    }
  } else {
    print(paste0(cluster, ": No markers to enrich"))
  }
}

# count number of DE genes after pseudobulk
upreg <- c()
downreg <- c()
for(i in unique(Idents(data))) {
  csv <- read.csv(file.path(dir, paste0("diff_bulk_", i, "_result.csv")))
  u <- csv |> filter(padj < 0.05 & log2FoldChange > 0.5) |> pull(gene) |> length()
  d <- csv |> filter(padj < 0.05 & log2FoldChange < -0.5) |> pull(gene) |> length()
  upreg <- c(upreg, u)
  downreg <- c(downreg, d)
}
de_n <- data.frame(
  cluster = unique(Idents(data)),
  upreg = upreg,
  downreg = downreg
)
write.csv(de_n, file.path(dir, "n_de_genes_pseudobulk_by_cluster.csv"))

# plot genes after pseudobulk
de_n <- read.csv(file.path(dir, "n_de_genes_pseudobulk_by_cluster.csv"))
de_n <- de_n |> mutate(total = upreg + downreg) |> arrange(total)
de_n$cluster = factor(de_n$cluster, levels = pull(de_n, cluster))
de_n <- de_n |> mutate(downreg = downreg * -1) |> 
  pivot_longer(cols=c(upreg, downreg), names_to = "type", values_to = "n") 
de_p <- de_n |>
  ggplot(aes(
    x = cluster,
    y = n,
    fill = type,
    label = abs(n)
  )) + theme_classic() +
  geom_col(data = subset(de_n, type == "upreg"), colour = "black") + coord_flip() + 
  geom_text(data = subset(de_n, type == "upreg"), nudge_y = -2) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) signif(abs(x), 3)) #+
#theme(panel.spacing.x = unit(0, "mm")) 
de_p <- de_p +geom_col(data = subset(de_n, type == "downreg"), colour = "black") + coord_flip() + 
  geom_text(data = subset(de_n, type == "downreg"), nudge_y = 2) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) signif(abs(x), 3), n.breaks = 14, limits = c(-70,70)) +
  scale_fill_manual(values=c("#00DAe0", "#E199FF"))
de_p
saveThePlot(path = file.path(dir, paste0("diff_bulk_n_barplot")), w = 6, h = 4)

# monocyte similarities
mono_shared <- list(
  Cluster_3 = read.csv(file.path(dir, paste0("diff_bulk_3_result.csv"))) |> filter(padj < 0.05) |> pull(gene),
  Cluster_6 = read.csv(file.path(dir, paste0("diff_bulk_6_result.csv"))) |> filter(padj < 0.05) |> pull(gene),
  Cluster_9 = read.csv(file.path(dir, paste0("diff_bulk_9_result.csv"))) |> filter(padj < 0.05) |> pull(gene)
)
ggvenn(
  mono_shared, 
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
saveThePlot(path = file.path(dir, paste0("venn_mono")), w = 8, h = 8)


###############################
## Augur
if(!require(Augur)) pak::pak("neurorestore/Augur")
library(Augur)

augur <- calculate_auc(data, label = "condition", cell_type_col = "seurat_clusters", n_threads = 8)
augur$AUC
write.csv(augur$AUC, file=file.path(dir, paste0("augur_auc.csv")))

important_genes <- augur$feature_importance |> arrange(cell_type, desc(importance)) |> filter(importance > 0.005)
write.csv(important_genes, file=file.path(dir, paste0("augur_important_genes.csv")))

# plot Augur
augur$AUC <- read.csv(file=file.path(dir, paste0("augur_auc.csv")))

augur$AUC |> mutate(cell_type = factor(cell_type, levels=unique(augur$AUC$cell_type))) |> ggplot(aes(
  x=cell_type,
  y=auc
)) + 
  geom_hline(yintercept=0.6, linetype='dashed', color='lightgray', size=0.5) +
  geom_segment(aes(x=cell_type, xend=cell_type, y=0.5, yend=auc)) +
  geom_point(color="orange", size=4) +
  geom_text(aes(label=signif(auc,2)), size = 2, position=position_nudge(y=0.02)) +
  scale_y_continuous(limits=c(0.5,0.8), expand=expansion(c(0,0))) +
  theme_classic()

saveThePlot(path = file.path(dir, paste0("Augur_Lollipop")), w = 5, h = 3)

##############################

## Differential abundance analysis with scProportionTest
if(!require(scProportionTest)) pak::pak("rpolicastro/scProportionTest")
library(scProportionTest)
prop_test <- sc_utils(data)
prop_test <- permutation_test(
  prop_test, cluster_identity = "seurat_clusters",
  sample_1 = "Control", sample_2 = "GBS",
  sample_identity = "condition"
)
write.csv(prop_test@results$permutation, file = file.path(dir, paste0("scProportionTest_Result_", Sys.Date(), ".csv")))

permutation_plot(prop_test)
saveThePlot(path = file.path(dir, paste0("scProportionTest")), w = 5, h = 4)

###################
#Per-sample PCA
# Create a matrix from our table of counts
pca_matrix <- AggregateExpression(data, group.by = "sample", 
                                  slot = "counts", 
                                  assays = "RNA",  
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000,
                                  return.seurat = F) 

# Perform the PCA
sample_pca <- prcomp(pca_matrix$RNA)

#extract eigenvalues
pc_eigenvalues <- sample_pca$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) |> 
  mutate(pct = variance/sum(variance)*100) |>
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

pc_eigenvalues |>
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# plot
pc_scores <- sample_pca$rotation |> as_tibble(rownames = "sample")
#pc_scores$condition <- c("GBS", "Control", "GBS", "Control", "Control", "GBS", "GBS")
pc_scores$condition <- c("GBS", "Control", "GBS", "Control", "Control", "GBS", "GBS")
pc_scores |> ggplot(aes(
  x = PC1, y = PC2,
  color = condition, label = sample)) +
  geom_point() + 
  geom_text(position = position_nudge(y = 0.05)) + 
  scale_color_manual(values=c("blue", "red")) +
  theme_bw()
saveThePlot(path = file.path(dir, paste0("SamplePCA")), w = 5, h = 4)

# top genes contributing to PC1 and PC2
pc_loadings <- sample_pca$x |> as_tibble(rownames = "gene")
top_genes <- pc_loadings |> 
  select(gene, PC1, PC2) |>
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") |>
  group_by(PC) |>
  arrange(desc(abs(loading))) |> 
  dplyr::slice(1:10) |>
  pull(gene) |> unique()

top_genes

top_loadings <- pc_loadings |> 
  filter(gene %in% top_genes)

ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "red") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            position = position_nudge(y = 30000), size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  theme_bw()
saveThePlot(path = file.path(dir, paste0("SamplePCA_LoadingsTop10")), w = 4.5, h = 4)

write.csv(tidy(sample_pca, matrix = "eigenvalues"), file=file.path(dir, "SamplePCA_Eigenvalues.csv")) # save PC variances
write.csv(tidy(sample_pca, matrix = "x"), file=file.path(dir, "SamplePCA_Loadings.csv")) # save variable loadings
write.csv(top_loadings, file=file.path(dir, "SamplePCA_Loadings_Top.csv")) # save variable loadings


#### sample correlations
bulk <- subset(data, seurat_clusters == 9)
bulk.subset <- AggregateExpression(bulk, group.by = c("condition", "sample"), 
                                   slot = "counts", 
                                   assays = "RNA",  
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000,
                                   return.seurat = F)
ccl4 <- bulk.subset$RNA["CCL4",]
ccl4 <- ccl4 |> as.data.frame() |> rownames_to_column(var = "sample")
ccl4 <- ccl4[4:7,] # GBS only
ccl4 <- ccl4 |> left_join(gbs_meta, by = join_by("sample" == "sample.condition"))

# stat
ccl4.cor <- cor(ccl4$ccl4, ccl4$csf.leukocytes, method = "pearson") # CSF leukocytes
cor(ccl4$ccl4, ccl4$csf.protein, method = "pearson") # CSF protein

# plot
ccl4 |> ggplot(aes(
  x = csf.leukocytes,
  y = ccl4
)) + geom_point(aes(color = severity)) + 
  geom_smooth(method = "lm") + 
  geom_text(aes(label = sample.id), position=position_nudge(y=150)) +
  theme_classic()
saveThePlot(path = file.path(dir, paste0("CCL4vsCSFLeuk")), w = 4, h = 3.5)

ccl4 |> ggplot(aes(
  x = csf.protein,
  y = ccl4
)) + geom_point(aes(color = severity)) + 
  #geom_smooth(method = "lm") + 
  geom_text(aes(label = sample.id), position=position_nudge(y=150)) +
  geom_smooth(method="lm", se = F, ) + 
  theme_classic()
saveThePlot(path = file.path(dir, paste0("CCL4vsCSFProt")), w = 4, h = 3.5)

# lollipop plots
ccl4 |> ggplot(aes(
  x=sample.id
)) + 
  geom_segment(aes(x=sample.id, xend=sample.id, y=0, yend=ccl4)) +
  geom_point(aes(y = ccl4, color = severity),  size=4) +
  #geom_point(aes(y = csf.leukocytes / coeff), color="tomato", size=4) +
  scale_color_manual(values=c("steelblue", "tomato")) +
  scale_y_continuous(expand=expansion(c(0,0.1))) +  #limits=c(0.5,0.8), 
  ylab("CCL4 (Normalized Expression)") +
  coord_flip() + scale_x_discrete(limits=rev) +
  theme_bw()
saveThePlot(path = file.path(dir, paste0("CCL4")), w = 5, h = 3)

ccl4 |> ggplot(aes(
  x=sample.id
)) + 
  geom_segment(aes(x=sample.id, xend=sample.id, y=0, yend=csf.leukocytes)) +
  geom_point(aes(y = csf.leukocytes, color = severity),  size=4) +
  #geom_point(aes(y = csf.leukocytes / coeff), color="tomato", size=4) +
  scale_color_manual(values=c("steelblue", "tomato")) +
  scale_y_continuous(expand=expansion(c(0,0.1))) +  #limits=c(0.5,0.8), 
  ylab("CSF Leukocytes (cells/cm^3)") +
  coord_flip() + scale_x_discrete(limits=rev) +
  theme_bw()
saveThePlot(path = file.path(dir, paste0("CSF_Leuk")), w = 5, h = 3)

###################
# comparison with publicly available viral infection data sets
# first, collect all GBS DE genes
upreg <- c()
downreg <- c()
for(i in unique(Idents(data))) {
  csv <- read.csv(file.path(dir, paste0("diff_bulk_", i, "_result.csv")))
  u <- csv |> filter(padj < 0.05 & log2FoldChange > 0.5) |> pull(gene)
  d <- csv |> filter(padj < 0.05 & log2FoldChange < -0.5) |> pull(gene)
  upreg <- c(upreg, u)
  downreg <- c(downreg, d)
}
bulk.de.cd8 <- read.csv(file.path(dir, "CD8_11_sig_differential_expression.csv")) |> 
  mutate(log2FoldChange = log2FoldChange * -1) |> # it's in the wrong order right now
  mutate(diff = case_when(
    padj < 0.05 & log2FoldChange > 0.5 ~ "upreg",
    padj < 0.05 & log2FoldChange < -0.5 ~ "downreg",
    padj >= 0.05 ~ "ns"
  ))
bulk.de.cd11b <- read.csv(file.path(dir, "CD11b_11_sig_differential_expression.csv")) |> 
  mutate(log2FoldChange = log2FoldChange * -1) |> # it's in the wrong order right now
  mutate(diff = case_when(
    padj < 0.05 & log2FoldChange > 0.5 ~ "upreg",
    padj < 0.05 & log2FoldChange < -0.5 ~ "downreg",
    padj >= 0.05 ~ "ns"
  ))
cd8.up <- bulk.de.cd8 |> filter(diff=="upreg") |> filter(!is.na(symbol) & symbol !="") |> pull(symbol)
cd8.down <- bulk.de.cd8 |> filter(diff=="downreg") |> filter(!is.na(symbol) & symbol !="") |> pull(symbol)
cd11b.up <- bulk.de.cd11b |> filter(diff=="upreg") |> filter(!is.na(symbol) & symbol !="") |> pull(symbol)
cd11b.down <- bulk.de.cd11b |> filter(diff=="downreg") |> filter(!is.na(symbol) & symbol !="") |> pull(symbol)

upreg <- c(upreg, cd8.up, cd11b.up) |> unique()
downreg <- c(downreg, cd8.up, cd11b.up) |> unique()

# get viral infection genes
viral.genes <- read.csv(file.path(dir, "viral_genes.csv"))|>  # upregulated genes common in COVID-19 and influenza
  separate_longer_delim(u.flu.covid, delim = ", ") |>
  separate_longer_delim(u.flu.only, delim = ", ") |>
  separate_longer_delim(u.covid.only, delim = ", ") |>
  separate_longer_delim(u.ebv, delim = ", ")

u.flu.covid <- viral.genes |> separate_longer_delim(u.flu.covid, delim = ", ") |> 
  filter(u.flu.covid != "") |> pull(u.flu.covid) |> unique()
u.flu.only <- viral.genes |> separate_longer_delim(u.flu.only, delim = ", ") |> 
  filter(u.flu.only != "") |> pull(u.flu.only) |> unique()
u.covid.only <- viral.genes |> separate_longer_delim(u.covid.only, delim = ", ") |> 
  filter(u.covid.only != "") |> pull(u.covid.only) |> unique()
u.ebv <- viral.genes |> separate_longer_delim(u.ebv, delim = ", ") |> 
  filter(u.ebv != "") |> pull(u.ebv) |> unique()
u.flu <- c(u.flu.only, u.flu.covid)
u.covid <- c(u.covid.only, u.flu.covid)

upset.genes <- list(GBS=upreg, Influenza = u.flu, `SARS-CoV-2` = u.covid, EBV = u.ebv)

# plot sets

if(!require(UpsetR)) install.packages("UpSetR")
library(UpSetR)

tiff(file.path(dir, "Viral_UpSet.tiff"), width=3000, height=3000, res=600)
upset(fromList(upset.genes),
      sets = c("EBV", "SARS-CoV-2", "Influenza", "GBS"),
      number.angles = 0,
      mb.ratio = c(0.67, 0.33),
      show.numbers=F,
      point.size = 3,
      line.size = 1,
      sets.bar.color = "steelblue", 
      main.bar.color = "grey",
      queries = list(
        list(query = intersects,
             params = list("GBS"),
             color = "#33b31a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "EBV", "SARS-CoV-2", "Influenza"), 
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "Influenza"), 
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "Influenza", "SARS-CoV-2"), 
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "SARS-CoV-2"), 
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "EBV"), 
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "EBV", "Influenza"),
             color = "#b3331a", active = TRUE),
        list(query = intersects,
             params = list("GBS", "SARS-CoV-2", "EBV"),
             color = "#b3331a", active = TRUE)
      ),
      #keep.order = T)
      order.by = "freq")
dev.off()

# output sets
gbs.viral.up <- intersect(upreg, c(u.flu, u.covid, u.ebv))
gbs.nonviral.up <- setdiff(upreg, c(u.flu, u.covid, u.ebv))

# enrichment analysis of GBS-specific genes
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023","KEGG_2016", "Panther_2016")

enriched.gbs <- enrichr(gbs.nonviral.up, dbs)
plotEnrich(enriched.gbs[[2]], showTerms = 10, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("Virus_GO_BP_Nonviral_upreg")), w = 8, h = 6)
write.csv(enriched.gbs[[2]], file = file.path(dir, paste0("Virus_GO_BP_Nonviral_upreg.csv")))
head(enriched.gbs[[2]], 10) |> ggplot(aes(x = Adjusted.P.value, y = Term)) +
  geom_dotplot()

# enrichment analysis of GBS genes also expressed in viral infection data sets
enriched.vir <- enrichr(gbs.viral.up, dbs)
plotEnrich(enriched.vir[[2]], showTerms = 10, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir, paste0("Virus_GO_BP_Viral_upreg")), w = 8, h = 6)
write.csv(enriched.vir[[2]], file = file.path(dir, paste0("Virus_GO_BP_Viral_upreg.csv")))
head(enriched.vir[[2]], 10) |> ggplot(aes(x = Adjusted.P.value, y = Term)) +
  geom_dotplot()

# plot pathway odds ratios against each other
# collect data
odds <- enriched.gbs[[2]] |> full_join(enriched.vir[[2]], by=join_by(Term))
odds.x <- odds |> filter(Adjusted.P.value.x < 0.05) |> #arrange(desc(Odds.Ratio.x)) |> head(20)# |>
  arrange(Adjusted.P.value.x) |> head(25)
odds.y <- odds |> filter(Adjusted.P.value.y < 0.05) |> filter(Term != "Interleukin-27-Mediated Signaling Pathway (GO:0070106)") |> # only one gene in this set ->  large outlier
  arrange(Adjusted.P.value.y) |> head(25)

# plot cleveland
odds.x |> ggplot() +
  geom_segment(aes(x=reorder(Term, Odds.Ratio.x), xend=Term, y=Odds.Ratio.x, yend=Odds.Ratio.y), color="grey") +
  geom_point(aes(x=Term, y=Odds.Ratio.x), color=rgb(0.2,0.7,0.1,0.5), size=2) + # GBS-specific in green
  geom_point(aes(x=Term, y=Odds.Ratio.y), color=rgb(0.7,0.2,0.1,0.5), size=2) + # GBS-viral in red
  coord_flip()+ 
  theme_classic() +
  xlab("") +
  ylab("Odds Ratio")
saveThePlot(path = file.path(dir, paste0("Viral_GBS-Specific_Cleveland")), w = 6, h = 4)

odds.y |> ggplot() +
  geom_segment(aes(x=reorder(Term, Odds.Ratio.y), xend=Term, y=Odds.Ratio.x, yend=Odds.Ratio.y), color="grey") +
  geom_point(aes(x=Term, y=Odds.Ratio.x), color=rgb(0.2,0.7,0.1,0.5), size=3) + # GBS-specific in green
  geom_point(aes(x=Term, y=Odds.Ratio.y), color=rgb(0.7,0.2,0.1,0.5), size=3) + # GBS-viral in red
  coord_flip()+ 
  theme_classic() +
  xlab("") +
  ylab("Odds Ratio")
saveThePlot(path = file.path(dir, paste0("Viral_GBS-Viral_Cleveland")), w = 6, h = 4)


###################
# DEBUG AND REPORTING
writeLines(capture.output(sessionInfo()), file.path(dir, paste0("sessionInfo_", Sys.Date(), ".txt")))
