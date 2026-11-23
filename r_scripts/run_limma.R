#!/usr/bin/env Rscript
# BioPrism - Differential Expression with Limma
# Usage: Rscript run_limma.R <expr_csv> <meta_csv> <output_csv> <case_group> <control_group>

suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Error: Insufficient arguments.\n")
  cat("Usage: Rscript run_limma.R <expr_csv> <meta_csv> <output_csv> <case_group> <control_group>\n")
  quit(status = 1)
}

expr_file <- args[1]
meta_file <- args[2]
out_file  <- args[3]
case_grp  <- args[4]
ctrl_grp  <- args[5]

cat(sprintf("BioPrism Limma Runner: Starting DE analysis\n"))
cat(sprintf("Expression file: %s\n", expr_file))
cat(sprintf("Metadata file: %s\n", meta_file))
cat(sprintf("Comparison: %s vs %s (Control)\n", case_grp, ctrl_grp))

tryCatch({
  # 1. Read Expression Matrix
  expr_data <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
  
  # 2. Read Metadata
  meta_data <- read.csv(meta_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(meta_data)[1:2] <- c("Sample", "Group")
  
  # Check matching samples
  common_samples <- intersect(colnames(expr_data), meta_data$Sample)
  if (length(common_samples) < 2) {
    stop("Less than 2 common samples found between expression data and metadata.")
  }
  
  # Align data
  meta_sub <- meta_data[meta_data$Sample %in% common_samples, ]
  meta_sub <- meta_sub[match(common_samples, meta_sub$Sample), ]
  expr_sub <- expr_data[, common_samples]
  
  # Filter for groups of interest
  valid_idx <- which(meta_sub$Group %in% c(case_grp, ctrl_grp))
  if (length(valid_idx) < 2) {
    stop("Insufficient samples belonging to the specified Case and Control groups.")
  }
  
  meta_sub <- meta_sub[valid_idx, ]
  expr_sub <- expr_sub[, valid_idx]
  
  # Count per group
  n_case <- sum(meta_sub$Group == case_grp)
  n_ctrl <- sum(meta_sub$Group == ctrl_grp)
  if (n_case < 2 || n_ctrl < 2) {
    stop(sprintf("Limma requires at least 2 replicates per group. Found Case: %d, Control: %d", n_case, n_ctrl))
  }
  
  # Ensure factors are syntactically valid R names for contrast
  group_levels <- c("Control", "Case")
  group_factor <- factor(ifelse(meta_sub$Group == case_grp, "Case", "Control"), levels = group_levels)
  
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- c("Control", "Case")
  
  # Fit linear model
  fit <- lmFit(as.matrix(expr_sub), design)
  
  # Compute contrast: Case - Control
  contrast_matrix <- makeContrasts(CaseVsControl = Case - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Extract all results
  res_table <- topTable(fit2, coef = "CaseVsControl", number = Inf, adjust.method = "BH", sort.by = "P")
  
  # Format output
  res_df <- data.frame(
    Gene = rownames(res_table),
    log2FoldChange = res_table$logFC,
    AveExpr = res_table$AveExpr,
    t_statistic = res_table$t,
    p_value = res_table$P.Value,
    adj_p_val = res_table$adj.P.Val,
    B_statistic = res_table$B,
    stringsAsFactors = FALSE
  )
  
  write.csv(res_df, out_file, row.names = FALSE)
  cat(sprintf("Successfully completed limma analysis. Processed %d genes. Results written to %s\n", nrow(res_df), out_file))
  quit(status = 0)
}, error = function(e) {
  cat(sprintf("Limma Execution Error: %s\n", e$message))
  quit(status = 2)
})
