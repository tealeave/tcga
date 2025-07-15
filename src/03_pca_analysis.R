# PCA Analysis Script
# Performs PCA analysis on raw counts and CPM-normalized data following blog methodology

suppressMessages({
  library(edgeR)
  library(dplyr)
})

cat("Loading PCA analysis functions...\n")

# Load preprocessed data
load(file.path(data_dir, "combined_counts.RData"))

# Function to select top variable genes
select_top_variable_genes <- function(count_matrix, n_genes = 1000) {
  cat("Selecting top", n_genes, "variable genes...\n")
  
  # Calculate gene variances
  gene_vars <- apply(count_matrix, 1, var)
  
  # Select top variable genes
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_genes])
  
  return(count_matrix[top_var_genes, ])
}

# Function to perform PCA analysis
perform_pca <- function(count_matrix, sample_metadata, analysis_type = "raw") {
  cat("Performing PCA analysis on", analysis_type, "data...\n")
  
  # Transpose matrix (samples as rows, genes as columns)
  pca_matrix <- t(count_matrix)
  
  # Perform PCA
  pca_result <- prcomp(pca_matrix, center = TRUE, scale. = FALSE)
  
  # Extract PC scores
  pc_scores <- data.frame(pca_result$x)
  pc_scores$sample_id <- rownames(pc_scores)
  pc_scores$cancer_type <- sample_metadata$cancer_type
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  # Store results
  pca_results <- list(
    pca_obj = pca_result,
    pc_scores = pc_scores,
    var_explained = var_explained,
    analysis_type = analysis_type
  )
  
  return(pca_results)
}

# Analysis 1: Raw counts PCA
cat("=== Raw Counts PCA Analysis ===\n")
raw_counts_subset <- select_top_variable_genes(combined_counts, 1000)
raw_pca <- perform_pca(raw_counts_subset, sample_metadata, "raw_counts")

# Analysis 2: CPM-normalized PCA
cat("=== CPM-Normalized PCA Analysis ===\n")
# Calculate CPM using edgeR
cpm_counts <- cpm(combined_counts, log = FALSE)
cpm_log_counts <- cpm(combined_counts, log = TRUE, prior.count = 1)

# Select top variable genes from CPM data
cpm_counts_subset <- select_top_variable_genes(cpm_log_counts, 1000)
cpm_pca <- perform_pca(cpm_counts_subset, sample_metadata, "cpm_normalized")

# Save PCA results
save(raw_pca, cpm_pca, file = file.path(data_dir, "pca_results.RData"))

# Save results tables
write.csv(raw_pca$pc_scores, file = file.path(tables_dir, "raw_pca_scores.csv"), row.names = FALSE)
write.csv(cpm_pca$pc_scores, file = file.path(tables_dir, "cpm_pca_scores.csv"), row.names = FALSE)

# Create variance explained summary
variance_summary <- data.frame(
  PC = paste0("PC", 1:10),
  raw_counts_var = raw_pca$var_explained[1:10],
  cpm_normalized_var = cpm_pca$var_explained[1:10]
)

write.csv(variance_summary, file = file.path(tables_dir, "variance_explained.csv"), row.names = FALSE)

cat("PCA analysis completed successfully!\n")
cat("Raw counts - PC1 variance explained:", round(raw_pca$var_explained[1], 2), "%\n")
cat("Raw counts - PC2 variance explained:", round(raw_pca$var_explained[2], 2), "%\n")
cat("CPM normalized - PC1 variance explained:", round(cpm_pca$var_explained[1], 2), "%\n")
cat("CPM normalized - PC2 variance explained:", round(cpm_pca$var_explained[2], 2), "%\n")