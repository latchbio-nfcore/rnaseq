#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

suppressPackageStartupMessages({
library(optparse)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
})

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=3       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

args <- commandArgs(trailingOnly = TRUE)

rds_path <- args[1]
condition_file <- args[2]
output_dir <- args[3]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the RDS file
se <- readRDS(rds_path)
counts <- assays(se)$counts

sampleTable <- read.csv(condition_file, row.names = 1)

# Filter out samples marked as "Ignore"
keep_samples <- rownames(sampleTable)[sampleTable$condition != "Ignore"]
sampleTable <- sampleTable[keep_samples, , drop = FALSE]
counts <- counts[, keep_samples]

print("Dimensions of counts:")
print(dim(counts))
print("Head of counts:")
print(head(counts))

print("Dimensions of sampleTable:")
print(dim(sampleTable))
print("Head of sampleTable:")
print(head(sampleTable))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = sampleTable, design = ~ condition)
dds <- DESeq(dds)

saveRDS(dds, file = file.path(output_dir, "dds.rds"))
res <- results(dds)
write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))

# Save size factor QC data to CSV
size_factors_data <- data.frame(Sample = colnames(counts(dds)), SizeFactor = sizeFactors(dds))
write.csv(size_factors_data, file = file.path(output_dir, "Size_Factor_QC.csv"))

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

# PCA data
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
write.csv(pcaData, file = file.path(output_dir, "PCA_data.csv"))

# Sample distance heatmap data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
write.csv(sampleDistMatrix, file = file.path(output_dir, "Sample_distance_matrix.csv"))

# Sample correlation data
vsd_cor <- cor(assay(vsd))
write.csv(vsd_cor, file = file.path(output_dir, "Sample_Correlation.csv"))

# Heatmap of the top 100 most expressed genes with Z scores
vsd_assay <- assay(vsd)
top_genes <- head(order(rowMeans(vsd_assay), decreasing = TRUE), 100)
heatmap_data <- vsd_assay[top_genes, ]

# Z-score heatmap for the top 100 most expressed genes
heatmap_data_zscores <- t(scale(t(heatmap_data), center = TRUE, scale = TRUE))
write.csv(heatmap_data_zscores, file = file.path(output_dir, "Counts_Heatmap_Zscores.csv"))
