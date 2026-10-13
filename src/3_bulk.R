library(tidyverse)
library(ggplot2)
library(enrichR)
library(ggvenn)

# load custom functions
source("src", "utils.R")

type <- "CD8" # CD4, CD8, CD11b
dir_ds <- file.path("DESeq", paste0("DESeq_results_", type))

# DESeq2
bulk.gbs <- read.csv(file.path(dir_ds, "10_results_annotated.csv")) # choose the annotated results file for DESeq2
bulk.gbs <- bulk.gbs |> mutate(log_padj = -log10(padj)) |> mutate(log_pvalue = -log10(pvalue)) |>
  mutate(log2FoldChange = log2FoldChange * -1) |> # input used GBS as reference group so has wrong direction - reverse the sign
  mutate(diff = case_when(
    padj < 0.05 & log2FoldChange > 0.5 ~ "upreg",
    padj < 0.05 & log2FoldChange < -0.5 ~ "downreg",
    padj >= 0.05 ~ "ns"
  ))

bulk.gbs$diff <- factor(bulk.gbs$diff, levels = c("upreg", "downreg", "ns"))
n_diff <- bulk.gbs |> group_by(diff) |> summarise(n=n())
write.csv(n_diff, file = file.path(dir_ds, "n.csv"))

# replace empty gene name with ensemble ID
bulk.gbs <- bulk.gbs |> filter(!is.na(symbol))
bulk.gbs$symbol[bulk.gbs$symbol==""] <- bulk.gbs$ensgene[bulk.gbs$symbol==""]

# regular volcano

bulk.gbs |> ggplot(aes(
  x = log2FoldChange,
  y = log_pvalue
)) + geom_point(aes(colour = diff, size = diff)) + 
  scale_colour_manual(values = c("red", "blue", "gray")) +
  scale_size_manual(values = c(2, 2, 1)) +
  geom_text(aes(label = symbol), data = bulk.gbs[bulk.gbs$padj < 0.05,], check_overlap = T, nudge_y = 0.5) +
  theme_classic()
saveThePlot(path = file.path(dir_ds, paste0("Volcano")), w = 8, h = 6)

# enhanced volcano
genes.to.plot <- c("KREMEN1", "BMAL2", "FAS", "MYO7B", "CREB5", "ADAMTS1", "LGR6", 
                   "EREG", "GAS6", "ADM", "IL15", "FN1", "ANPEP", "SDC4", "IFIT3", 
                   "FCGBP", "NEIL1", "MYOM2", "BAIAP3", "ASAP2", "SERPING1",
                   "IFNB1", "IFIH1", "IFIT1", "IFIT2", "IFITM3", "GBP2", "OAS1", "OAS2",
                   "CXCL9", "CXCL10", "IL6", "IL31RA", "MARCO", "CCR5", "CLU", "STAT1") # CD11b

EnhancedVolcano(bulk.gbs,
                lab = bulk.gbs$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0, max(-log10(bulk.gbs[['pvalue']]), na.rm = TRUE) + 0.5),
                xlim = c(min(bulk.gbs[['log2FoldChange']], na.rm = TRUE) - 0.5, 
                         max(bulk.gbs[['log2FoldChange']], na.rm = TRUE) + 0.5),
                title = paste0(type, "+ cells"),
                subtitle = "",
                drawConnectors = T,
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                arrowheads = F,
                selectLab = genes.to.plot,
                max.overlaps = 30)
saveThePlot(path = file.path(dir_ds, paste0("Volcano_enhanced_", Sys.Date())), w = 9, h = 9)


upreg <- bulk.gbs |> filter(diff == "upreg") |> pull(symbol) |> na.omit()
upreg <- upreg[upreg!= ""]
downreg <- bulk.gbs |> filter(diff == "downreg") |> pull(symbol) |> na.omit()
downreg <- downreg[downreg!= ""]
upreg_ens <- bulk.gbs |> filter(diff == "upreg") |> pull(ensgene) |> na.omit()
upreg_ens <- upreg_ens[upreg_ens!= ""]
downreg_ens <- bulk.gbs |> filter(diff == "downreg") |> pull(ensgene) |> na.omit()
downreg_ens <- downreg_ens[downreg_ens!= ""]

# heatmap across individual samples
groups <- read.csv(file.path("Bulk RNA Seq", "JoeSamples", "misc", paste0("metadata_stratified_", type, ".csv")))
hm.data <- read.csv(file.path(dir_trap, "1_unfiltered_tpm_data.csv")) |>#"3_donor_wise_tpm.csv")) |>
  filter(Gene.ID %in% c(upreg_ens, downreg_ens)) |>
  pivot_longer(cols = starts_with("X356"), names_to = "sample_id", values_to = "counts") |>
  mutate(sample_id = gsub("X", "", sample_id)) |> mutate(sample_id = gsub("\\.", "-", sample_id)) |>
  left_join(groups, by = join_by(sample_id == Sample)) |>
  mutate(Group = case_when(
    Group == "HC" ~ "Control",
    Group == "GBS" ~ "GBS"
  )) |>
  group_by(Gene.Name) |> mutate(counts = log2(counts/mean(counts))) |> ungroup() |>
  filter(!is.na(Gene.Name)) |>
  left_join(bulk.gbs, by=join_by(Gene.ID==ensgene))
hm.data <- hm.data |> mutate(symbol = case_when(
  is.na(symbol) ~ Gene.ID,
  !is.na(symbol) ~ symbol
))
hm.data$counts[is.infinite(hm.data$counts)] <- NA
hm.data$sample_id <- factor(hm.data$sample_id, levels=unique(hm.data$sample_id[order(hm.data$Group, hm.data$sample_id)]), ordered=T)
hm.data$symbol <- factor(hm.data$symbol, levels=unique(hm.data$symbol[order(hm.data$diff, hm.data$symbol)]), ordered=T)

hm.data |> ggplot(aes(
  x = sample_id,
  y = symbol,
  fill = counts
)) + geom_tile() + 
  scale_fill_distiller(palette = "RdBu", type="div", limit = max(abs(hm.data$counts), na.rm=T) * c(-1, 1)) + 
  scale_x_discrete(position = "top") + scale_y_discrete(limits=rev) +
  facet_wrap(~Group, scales = "free_x") +
  labs(fill = "Normalized log2(TPM)") +
  theme_classic() +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
if(type == "CD11b") h <- 16 else h <- 7
saveThePlot(path = file.path(dir_ds, paste0("PerSample_Heatmap")), w = 3.5, h = h)
hm.data |> group_by(diff, Gene.Name) |> summarise() |> summarise(n=n())

# enrichR analysis
# enrich <- read.csv(file.choose())
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023","KEGG_2016", "Panther_2016")
#enriched <- enrichr(as.vector(enrich$symbol), dbs) # to take top 10 genes per cluster regardless of significance

#both upreg and downreg
enriched <- enrichr(c(upreg, downreg), dbs)
enriched[["GO_Molecular_Function_2023"]] #<- enriched[["GO_Molecular_Function_2023"]] |> filter(Adjusted.P.value < 0.05)
enriched[["GO_Biological_Process_2023"]] #<- enriched[["GO_Biological_Process_2023"]] |> filter(Adjusted.P.value < 0.05)
enriched[["KEGG_2016"]] #<- enriched[["KEGG_2016"]] |>  filter(Adjusted.P.value < 0.05)
enriched[["Panther_2016"]] #<- enriched[["Panther_2016"]] |> filter(Adjusted.P.value < 0.05)
write.csv(enriched[["GO_Molecular_Function_2023"]], file = file.path(dir_ds, "GO_MF_reg.csv"))
write.csv(enriched[["GO_Biological_Process_2023"]], file = file.path(dir_ds, "GO_BP_reg.csv"))
write.csv(enriched[["KEGG_2016"]], file = file.path(dir_ds, "KEGG_reg.csv"))
write.csv(enriched[["Panther_2016"]], file = file.path(dir_ds, "Panther_reg.csv"))
plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_MF_reg")), w = 8, h = 6)
plotEnrich(enriched[[2]], showTerms = 9, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_BP_reg")), w = 9, h = 2.8)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
plotEnrich(enriched[[4]] |> filter(P.value < 0.05), 
           showTerms = 10, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("Panther_reg")), w = 9, h = 2)

#upreg
enriched <- enrichr(upreg, dbs)#enrichr(c(upreg, downreg), dbs)
enriched[["GO_Molecular_Function_2023"]]
enriched[["GO_Biological_Process_2023"]]
enriched[["KEGG_2016"]]
enriched[["Panther_2016"]]
plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_MF_upreg")), w = 8, h = 6)
plotEnrich(enriched[[2]], showTerms = 10, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_BP_upreg")), w = 8, h = 6)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
plotEnrich(enriched[[4]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")

#downreg
enriched <- enrichr(downreg, dbs)#enrichr(c(upreg, downreg), dbs)
enriched[["GO_Molecular_Function_2023"]]
enriched[["GO_Biological_Process_2023"]]
enriched[["KEGG_2016"]]
enriched[["Panther_2016"]]
plotEnrich(enriched[[1]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_MF_downreg")), w = 8, h = 6)
plotEnrich(enriched[[2]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
saveThePlot(path = file.path(dir_ds, paste0("GO_BP_downreg")), w = 8, h = 6)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")
plotEnrich(enriched[[4]], showTerms = 20, numChar = 60, y = "Ratio", orderBy = "P.value")



# venn diagram of similarities between compartments
dir <- file.path("Bulk RNA Seq", "JoeSamples", "DESeq")
#upreg
bulk_shared <- list(
  CD4 = read.csv(file.path(dir, "DESeq_results_CD4", "11_sig_differential_expression.csv")) |> 
    mutate(log2FoldChange = log2FoldChange * -1) |>
    filter(log2FoldChange > 0.5) |> pull(symbol),
  CD8 = read.csv(file.path(dir, "DESeq_results_CD8", "11_sig_differential_expression.csv")) |> 
    mutate(log2FoldChange = log2FoldChange * -1) |>
    filter(log2FoldChange > 0.5) |> pull(symbol),
  CD11b = read.csv(file.path(dir, "DESeq_results_CD11b", "11_sig_differential_expression.csv")) |>
    mutate(log2FoldChange = log2FoldChange * -1) |>
    filter(log2FoldChange > 0.5) |> pull(symbol)
)
ggvenn(
  bulk_shared, 
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4#,
  #auto_scale = TRUE
)
saveThePlot(path = file.path(dir, paste0("venn_bulk_upreg")), w = 4, h = 4)

#downreg
bulk_shared <- list(
  CD4 = read.csv(file.path(dir, "DESeq_results_CD4", "11_sig_differential_expression.csv")) |> 
   mutate(log2FoldChange = log2FoldChange * -1) |>
   filter(log2FoldChange < -0.5) |> pull(symbol),
  CD8 = read.csv(file.path(dir, "DESeq_results_CD8", "11_sig_differential_expression.csv")) |> 
    mutate(log2FoldChange = log2FoldChange * -1) |>
    filter(log2FoldChange < -0.5) |> pull(symbol),
  CD11b = read.csv(file.path(dir, "DESeq_results_CD11b", "11_sig_differential_expression.csv")) |>
    mutate(log2FoldChange = log2FoldChange * -1) |>
    filter(log2FoldChange < -0.5) |> pull(symbol)
)
ggvenn(
  bulk_shared, 
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = TRUE
)
saveThePlot(path = file.path(dir, paste0("venn_bulk_downreg")), w = 4, h = 4)
