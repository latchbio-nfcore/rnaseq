#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)
file_list_path <- args[1]
tx2gene_path <- args[2]
output_dir <- args[3]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the list of rsem files
file_list <- readLines(file_list_path)
print("RSEM files:")
print(file_list)

# Extract sample names from file paths
sample_names <- sub("\\.isoforms\\.results$", "", basename(file_list))

# Read the tx2gene file
tx2gene <- read.table(tx2gene_path, header = FALSE, stringsAsFactors = FALSE)
tx2gene <- tx2gene[, 1:2]

# Import the RSEM data using tximport with lengthScaledTPM
txi.rsem <- tximport(file_list, type = "rsem",
                     txIn = TRUE, txOut = FALSE,
                     tx2gene = tx2gene,
                     countsFromAbundance = "lengthScaledTPM")
txi.rsem$length[txi.rsem$length <= 0] <- 1

# Prepare the counts data frame with proper row and column names
counts_df <- as.data.frame(txi.rsem$counts)
colnames(counts_df) <- sample_names
counts_df$gene_id <- rownames(counts_df)
counts_df <- counts_df[, c("gene_id", sample_names)]

scaled_counts_file <- file.path(output_dir, "deseq2_rsem_length_scaled_counts.tsv")
write.csv(counts_df, file = scaled_counts_file, quote = FALSE, row.names = FALSE, sep = "\t")

print(paste("Length-scaled TPM counts have been written to:", scaled_counts_file))
