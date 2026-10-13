library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)

# load custom functions
source("src", "utils.R")

drg <- read_excel(file.path("data", "human_DRG_sn_V1.xlsx"))
dir <- file.path("outs", "Analysis_20241007", "Interactome")

# get DE genes to use
cluster <- 7
de.genes <- read.csv(file.path(dir, paste0("diff_bulk_", cluster, "_result.csv"))) |> 
  filter(padj < 0.05 & log2FoldChange > 0) |> pull(gene)

# extract expression info
lig <- AggregateExpression(
  data,
  assays = "RNA",
  features = de.genes,
  return.seurat = FALSE,
  group.by = c("seurat_clusters", "condition"),
  add.ident = NULL,
  #layer = "data",
  verbose = TRUE
) 

lig <- lig$RNA |> as.data.frame() |> select(starts_with(paste0("g", cluster)))
write.csv(lig, file=file.path(dir, cluster, "_lig_out.csv"))

# load receptor data from single-nuc
scdata <- readRDS(file = file.path("data", "singleNuc_hDRG.rds"))

msch <- subset(scdata, idents = "Mschwann") # limit to myelinating schwann
nmsch <- subset(scdata, idents = "NMSchwann") # limit to non-myelinating schwann

rec <- AggregateExpression(
  scdata,
  assays = "RNA",
  #features = de.genes,
  return.seurat = FALSE,
  group.by = "clusLabs",
  add.ident = NULL,
  #layer = "data",
  verbose = TRUE
) 
rec <- rec$RNA |> as.data.frame()
write.csv(rec, file=file.path(dir, "rec_out.csv"))

#### GBS: ALLUVIAL PLOTS ###########
library(ggalluvial)
library(readxl)
cell.toward <- "Neurons" # change as needed
colname <- sym(paste0("Rec_", cell.toward))
dir.i <- file.path("Analysis_20241007", "Interactome", "Cluster 6 Intermediate monocyte", cell.toward) # change cluster as needed
interactome_filter <- read_excel(file.path(dir.i, "Ligand_to_Receptor_RankedData(1).xlsx")) |>
#interactome_filter <- read_excel(file.path(dir.i, "Ligand_to_Receptor_RankedData_RecCatLimit.xlsx")) |>
  select(Ligand, lig_final_category, Receptor, rec_final_category, !!colname)

# load category colors
cat.color <- read_excel(file.path("data", "Category and Hex Color codes.xlsx")) |>
  mutate(Types = as.factor(Types))
#alluvial plot
int.long <- interactome_filter |> 
  mutate(Interaction = c(1:nrow(interactome_filter))) |> 
  pivot_longer(cols=c(Ligand, Receptor), names_to = "Molecule_Type", values_to = "Gene") |>
  mutate(Category = case_when(
    Molecule_Type == "Ligand" ~ lig_final_category,
    Molecule_Type == "Receptor" ~ rec_final_category
  )) |>
  select(-lig_final_category, -rec_final_category)

int.long |> 
  ggplot(aes(x = Molecule_Type, stratum = Gene, alluvium = Interaction, y = !!colname)) +
  geom_alluvium(aes(fill = Category), show.legend = FALSE) +
  geom_stratum(aes(fill = Category), show.legend = FALSE) +
  scale_fill_discrete(drop = T, limits = levels(cat.color$Types), 
                      type = setNames(cat.color$hexcodes, cat.color$Types)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme_void()

saveThePlot(path = file.path(dir.i, paste0("alluvial_SC_", "CD11b", "_", cell.toward, "_nTPM")), w = 4, h = 6) # 5 x 9 for CD11b, 4 x 3 for cluster 9

# GSEA OF RECEPTORS
markers_dea <- interactome_filter
if(dim(markers_dea)[1] > 0) { # check that df is not empty
  dbs <- listEnrichrDbs()
  dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "KEGG_2016", "Panther_2016")
  enriched <- enrichr(markers_dea$Receptor, dbs)
  enriched[["GO_Molecular_Function_2023"]]
  enriched[["GO_Biological_Process_2023"]]
  enriched[["KEGG_2016"]]
  enriched[["Panther_2016"]]
  if(dim(enriched[["GO_Molecular_Function_2023"]])[1] > 0) {
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "Combined.Score")
    saveThePlot(path = file.path(dir.i, "interactome_GO_MF_"), w = 8, h = 6)
    write.csv(enriched[[1]], file = file.path(dir.i, "interactome_GO_MF.csv"))
  }
  if(dim(enriched[["GO_Biological_Process_2023"]])[1] > 0) {
    plotEnrich(enriched[[2]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "Combined.Score")
    saveThePlot(path = file.path(dir.i,"interactome_GO_BP_"), w = 8, h = 6)
    write.csv(enriched[[2]], file = file.path(dir.i, paste0("diff_bulk_", cluster, "_GO_BP.csv")))
  }
  if(dim(enriched[["KEGG_2016"]])[1] > 0) {
    plotEnrich(enriched[[3]], showTerms = 50, numChar = 60, y = "Ratio", orderBy = "Combined.Score")
    saveThePlot(path = file.path(dir.i, "interactome_KEGG_"), w = 8, h = 6)
    write.csv(enriched[[3]], file = file.path(dir.i, "interactome_KEGG.csv"))
  }
  if(dim(enriched[["Panther_2016"]])[1] > 0) {
    plotEnrich(enriched[[4]], showTerms = 5, numChar = 60, y = "Ratio", orderBy = "Combined.Score")
    saveThePlot(path = file.path(dir.i, "interactome_Panther_"), w = 8, h = 6)
    write.csv(enriched[[4]], file = file.path(dir.i, "interactome_Panther.csv"))
  }
} else {
  print(paste0(cluster, ": No markers to enrich"))
}


# limit interactions to specific categories of interest
library(readxl)
dir <- file.path("Analysis_20241007", "Interactome", "Cluster 13 Non-classical monocytes", "MSchwann")
ranked <- read_xlsx(file.path(dir, "Ligand_to_Receptor_RankedData(1).xlsx"))
ranked |> filter(
  lig_final_category != "extracelluar matrix",
  rec_final_category %in% c("defense/immunity",
                            "G-protein coupled receptor",
                            "metabolite interconversion enzyme",
                            "ion channel",
                            "other transmembrane signal receptor",
                            "serine/threonine protein kinase receptor",
                            "tyrosine protein kinase receptor")
)

# count categories
ranked |> group_by(rec_final_category) |> summarize(n=n()) |> print(n=50)

# summary
dir <- file.path("Analysis_20241007", "Interactome")

int_sum <- read_xlsx(file.path(dir, "GBS_Interactome_Summary.xlsx"))
POP_REC <- "Neuron" # Mschwann or Neuron

#heatmap
int_sum <- int_sum |> filter(Population_Rec == POP_REC) |>
  select(-`Limited to categories`, -Population_Rec) |>
  column_to_rownames("Population_Lig") |>
  filter(if_any(!matches("ID"), ~ . != 0)) |>
  #select_if(function(col) mean(col) > 0) |> # if loading fresh, use this for Mschwann to get all the columns/categories with at least one gene across any population.
  select(all_of(c(cats_to_use))) |> # for neurons only after running Mschwann
  as.matrix()
int_sum

#cats_to_use <- colnames(int_sum) # Mschwann only- set the categories that Mschwann used in case you run Neuron afterwards

library("pheatmap")
hm <- pheatmap(int_sum, cluster_rows=F, show_rownames=TRUE,
               border_color = "black", number_color,
               cluster_cols=F, #scale = "row",
               gaps_row = 1, gaps_col = c(1,7),
               cellwidth = 12, cellheight = 12,
               angle_col = "45",
               display_numbers = T, number_format = "%i",
               filename = file.path(dir, paste0(
                 "interactome_heatmap_", POP_REC, "_", Sys.Date(), ".pdf")
               ))#, 
#annotation_col=df, 
#annotation_row=df2)
hm
save_pheatmap_pdf(hm, file.path(dir, paste0("interactome_heatmap_", POP_REC, "_", Sys.Date(), ".pdf")),  
                  width = 10, height = 6)