#define custom functions
saveThePlot <- function(path, plot = last_plot(), w = 8, h = 6, u = "in") {
  for(fileformat in c("pdf", "tiff")) {
    ggsave(plot = plot, filename = paste0(path, "_", Sys.Date(), ".", fileformat), 
           device = fileformat, width = w, height = h, units = u)   
  }
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

DESeq2DETest <- function(
    data.use,
    cells.1,
    cells.2,
    verbose = TRUE,
    ...
) {
  if (!PackageCheck('DESeq2', error = FALSE)) {
    stop("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
  }
  CheckDots(..., fxns = 'DESeq2::results')
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05,
    ...
  )
  to.return <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  return(to.return)
}