#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(plotly)
  library(heatmaply)
  library(EnhancedVolcano)
})

args <- commandArgs(trailingOnly = TRUE)

rds_path <- args[1]
condition_file <- args[2]
output_dir <- args[3]

files_dir <- file.path(output_dir, "Files")
graphs_dir <- file.path(output_dir, "Graphs")
contrast_files_dir <- file.path(files_dir, "Contrast")
contrast_graphs_dir <- file.path(graphs_dir, "Contrast")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(files_dir)) {
  dir.create(files_dir, recursive = TRUE)
}
if (!dir.exists(graphs_dir)) {
  dir.create(graphs_dir, recursive = TRUE)
}
if (!dir.exists(contrast_files_dir)) {
  dir.create(contrast_files_dir, recursive = TRUE)
}
if (!dir.exists(contrast_graphs_dir)) {
  dir.create(contrast_graphs_dir, recursive = TRUE)
}

# Load the RDS file
se <- readRDS(rds_path)
counts <- assays(se)$counts

sampleTable <- read.csv(condition_file, row.names = 1)

# Filter out samples marked as "Ignore"
keep_samples <- rownames(sampleTable)[sampleTable$condition != "Ignore"]
sampleTable <- sampleTable[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]

sampleTable$condition <- trimws(tolower(sampleTable$condition))

# Set the reference level
conditions <- unique(sampleTable$condition)

if ("control" %in% conditions) {
  # If 'control' exists, make it the reference
  sampleTable$condition <- relevel(factor(sampleTable$condition), ref = "control")
} else {
  # Otherwise, use the first alphabetical condition as reference
  sampleTable$condition <- factor(sampleTable$condition, levels = sort(conditions))
}

print("Dimensions of counts:")
print(dim(counts))
print("Head of counts:")
print(head(counts))

print("Dimensions of sampleTable:")
print(dim(sampleTable))
print("Head of sampleTable:")
print(head(sampleTable))

print("Levels of condition factor:")
print(levels(sampleTable$condition))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sampleTable,
                              design = ~ condition)
dds <- DESeq(dds)

saveRDS(dds, file = file.path(output_dir, "dds.rds"))
res <- results(dds)
write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))

# Save size factor QC data to CSV
size_factors_data <- data.frame(Sample = colnames(counts(dds)), SizeFactor = sizeFactors(dds))
write.csv(size_factors_data, file = file.path(files_dir, "Size_Factor_QC.csv"))

# Perform VST transformation for PCA and heatmap
vsd <- tryCatch(
  {
    vst(dds)
  },
  error = function(err) {
    message("Failed naive vst transform, attempting raw varianceStabilizingTransformation")
    varianceStabilizingTransformation(dds)
  }
)
saveRDS(vsd, file = file.path(output_dir, "vsd.rds"))

# PCA plot
print("Creating PCA plots")
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png(filename = file.path(graphs_dir, "PCA_plot.png"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
dev.off()
write.csv(pcaData, file = file.path(files_dir, "PCA_data.csv"))

pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, text = rownames(pcaData))) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggplotly(pca_plot) %>%
  htmlwidgets::saveWidget(file.path(graphs_dir, "PCA_plot.html"))

# Sample distance heatmap
print("Creating sample distance heatmaps")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
write.csv(sampleDistMatrix, file = file.path(files_dir, "Sample_distance_matrix.csv"))
png(filename = file.path(graphs_dir, "Sample_distance_heatmap.png"))
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)
dev.off()

heatmaply(sampleDistMatrix) %>%
  htmlwidgets::saveWidget(file.path(graphs_dir, "Sample_distance_heatmap.html"))

# Sample correlation heatmap
print("Creating sample correlation heatmaps")
vsd_cor <- cor(assay(vsd))
write.csv(vsd_cor, file = file.path(files_dir, "Sample_Correlation.csv"))
png(filename = file.path(graphs_dir, "Sample_correlation_heatmap.png"))
pheatmap(vsd_cor)
dev.off()

heatmaply(vsd_cor) %>%
  htmlwidgets::saveWidget(file.path(graphs_dir, "Sample_correlation_heatmap.html"))

# Heatmap of the top 100 most expressed genes with Z scores
print("Creating sample Z-score heatmap for top 100 most expressed genes")
vsd_assay <- assay(vsd)
top_genes <- head(order(rowMeans(vsd_assay), decreasing = TRUE), 100)
heatmap_data <- vsd_assay[top_genes, ]

heatmap_data_zscores <- t(scale(t(heatmap_data), center = TRUE, scale = TRUE))
write.csv(heatmap_data_zscores, file = file.path(files_dir, "Counts_Heatmap_Zscores.csv"))
png(filename = file.path(graphs_dir, "Counts_Heatmap_Zscores.png"))
pheatmap(heatmap_data_zscores)
dev.off()

heatmaply(heatmap_data_zscores) %>%
  htmlwidgets::saveWidget(file.path(graphs_dir, "Counts_Heatmap_Zscores.html"))

# Contrast compute for all conditions
print("Creating contrast files for all conditions")
sampleTable[["condition"]] <- as.factor(sampleTable[["condition"]])
ls <- levels(sampleTable[["condition"]])

if (length(ls) < 2) {
  stop(sprintf("Not enough levels in column %s for comparison", "condition"))
}

for (i in 1:length(ls)) {
  l1 <- ls[[i]]
  for (j in 1:length(ls)) {
    l2 <- ls[[j]]
    if (l1 == l2) {
      next
    }

    tryCatch({
      full <- sprintf("%s_vs_%s_%s", l1, l2, "condition")
      message(sprintf("Generating MA, and Volcano Plot for %s vs %s", l1, l2))

      res <- results(dds, contrast = c("condition", l1, l2))

      write.csv(as.data.frame(res), file = file.path(contrast_files_dir, sprintf("%s.csv", full)))
      dir.create(file.path(contrast_graphs_dir, full), recursive = TRUE)

      lfc <- lfcShrink(dds, res = res, type = "ashr")

      png(file = file.path(contrast_graphs_dir, full, "MA.png"), width = 960, height = 540)
      ma_plot <- plotMA(lfc, ylim = c(-5, 5), main = paste(l1, l2, sep = " vs "))
      print(ma_plot)
      dev.off()

      voc2 <- EnhancedVolcano(
        lfc,
        lab = rownames(lfc),
        drawConnectors = TRUE,
        x = "log2FoldChange",
        y = "padj",
        title = sprintf("%s vs %s", l1, l2),
        subtitle = "",
        legendPosition = "none",
        widthConnectors = 0.5
      )
      png(file = file.path(contrast_graphs_dir, full, "Volcano.png"), width = 960, height = 540)
      print(voc2)
      dev.off()

    }, error = function(err) {
      message(sprintf("Failed to generate plots for %s vs %s: %s", l1, l2, err$message))
    })
  }
}
