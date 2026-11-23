#!/usr/bin/env Rscript
# BioPrism - Differential Expression with DESeq2 for Raw Count Data
# Usage: Rscript run_deseq2.R <counts_csv> <meta_csv> <output_csv> <case_group> <control_group>

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Error: Insufficient arguments.\n")
  cat("Usage: Rscript run_deseq2.R <counts_csv> <meta_csv> <output_csv> <case_group> <control_group>\n")
  quit(status = 1)
}

counts_file <- args[1]
meta_file   <- args[2]
out_file    <- args[3]
case_grp    <- args[4]
ctrl_grp    <- args[5]

cat(sprintf("BioPrism DESeq2 Runner: Starting RNA-seq count DE analysis\n"))
cat(sprintf("Counts file: %s\n", counts_file))
cat(sprintf("Metadata file: %s\n", meta_file))
cat(sprintf("Comparison: %s vs %s (Control)\n", case_grp, ctrl_grp))

tryCatch({
  # 1. Read Counts Matrix (round to integer if non-integer)
  counts_data <- read.csv(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
  counts_matrix <- round(as.matrix(counts_data))
  storage.mode(counts_matrix) <- "integer"
  
  # 2. Read Metadata
  meta_data <- read.csv(meta_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(meta_data)[1:2] <- c("Sample", "Group")
  
  # Align samples
  common_samples <- intersect(colnames(counts_matrix), meta_data$Sample)
  if (length(common_samples) < 2) {
    stop("Less than 2 common samples between count matrix and metadata.")
  }
  
  meta_sub <- meta_data[meta_data$Sample %in% common_samples, ]
  meta_sub <- meta_sub[match(common_samples, meta_sub$Sample), ]
  counts_sub <- counts_matrix[, common_samples]
  
  # Filter groups
  valid_idx <- which(meta_sub$Group %in% c(case_grp, ctrl_grp))
  if (length(valid_idx) < 2) {
    stop("Insufficient samples for specified Case and Control groups.")
  }
  
  meta_sub <- meta_sub[valid_idx, ]
  counts_sub <- counts_sub[, valid_idx]
  
  n_case <- sum(meta_sub$Group == case_grp)
  n_ctrl <- sum(meta_sub$Group == ctrl_grp)
  if (n_case < 2 || n_ctrl < 2) {
    stop(sprintf("DESeq2 requires at least 2 replicates per condition. Found Case: %d, Control: %d", n_case, n_ctrl))
  }
  
  meta_sub$condition <- factor(ifelse(meta_sub$Group == case_grp, "Case", "Control"), levels = c("Control", "Case"))
  rownames(meta_sub) <- meta_sub$Sample
  
  # Filter out low count genes (at least 10 counts across all samples)
  keep <- rowSums(counts_sub) >= 10
  counts_sub <- counts_sub[keep, ]
  
  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = ~ condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Results contrast: condition Case vs Control
  res <- results(dds, contrast = c("condition", "Case", "Control"))
  
  # Format output
  res_df <- data.frame(
    Gene = rownames(res),
    log2FoldChange = res$log2FoldChange,
    AveExpr = log2(res$baseMean + 1),
    stat = res$stat,
    p_value = res$pvalue,
    adj_p_val = res$padj,
    baseMean = res$baseMean,
    stringsAsFactors = FALSE
  )
  
  # Handle NA adj_p_val (DESeq2 sets NA for low-count or extreme outliers)
  res_df$p_value[is.na(res_df$p_value)] <- 1.0
  res_df$adj_p_val[is.na(res_df$adj_p_val)] <- 1.0
  
  # Sort by p_value
  res_df <- res_df[order(res_df$p_value), ]
  
  write.csv(res_df, out_file, row.names = FALSE)
  cat(sprintf("Successfully completed DESeq2 analysis. Processed %d genes. Results written to %s\n", nrow(res_df), out_file))
  quit(status = 0)
}, error = function(e) {
  cat(sprintf("DESeq2 Execution Error: %s\n", e$message))
  quit(status = 2)
})
