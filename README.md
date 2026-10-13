# gbs-pbmc-rnaseq

This repository archives code used to generate the downstream data analysis and most figures related to the research project "Transcriptomic signatures of early, untreated Guillain-Barré syndrome: a preliminary study of predicted immune cell-peripheral nerve interactions."

The manuscript is preprinted at bioRxiv: <https://www.biorxiv.org/content/10.1101/2025.09.24.678413>

# Code

The following code files are included in src/:

- **1_sc_object_creation.R**: Code used to take .h5 CellBender outputs, perform initial QC and processing steps on a per-sample basis, and integrate samples into a single Seurat object for downstream analysis.

- **2_sc_downstream.R:** Takes the output of (1) for downstream cell type annotations, differential abundance and expression analysis, and other visualisations and analysis of the single-cell data referenced in the paper.

- **3_bulk.R:** Code for bulk RNA sequencing differential expression analysis. Can be run independently.

- **4_interactome.R:** Code for ligand-receptor interactome analysis. Takes inputs from both (2) and (3).

- **CCL4_Ca.Rmd:** R markdown workbook for analysis of calcium imaging data presented in Figure 7c onwards.

- **CaImagingTraces.Rmd:** Code for plotting average calcium imaging traces from raw data as shown in Figure 7b.

- **utils.R:** Custom functions, mainly for ease of saving plot outputs from ggplot or pheatmap.

# Data

- Single-cell: The processed Seurat .rds file that can be used as input in (2) is available at GEO (GSE304872).

- Bulk-seq: Normalized counts for CD4+, CD8+, and CD11b+ cells used as input for (3) are available at GEO (GSE304871).

# Acknowledgements

This work is the result of a collaboration between UT Dallas and UAB. Additional code contributions were made by Ishwarya Sankaranarayanan and Katherin Gabriel.
