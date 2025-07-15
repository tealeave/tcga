# Marker Gene Analysis Script
# Analyzes marker genes (NAPSA, TFF1, TP63, KRT5) following blog methodology

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(edgeR)
  library(tidyr)
})

cat("Loading marker gene analysis functions...\n")

# Load data
load(file.path(data_dir, "combined_counts.RData"))
load(file.path(data_dir, "pca_results.RData"))

# Define marker genes from the blog
marker_genes <- c("NAPSA", "TFF1", "TP63", "KRT5")

# Function to analyze marker gene expression
analyze_marker_genes <- function(count_matrix, sample_metadata, marker_genes) {
  cat("Analyzing marker genes:", paste(marker_genes, collapse = ", "), "\n")
  
  # Calculate CPM
  cpm_counts <- cpm(count_matrix, log = TRUE, prior.count = 1)
  
  # Check which marker genes are present
  present_markers <- intersect(marker_genes, rownames(cpm_counts))
  missing_markers <- setdiff(marker_genes, rownames(cpm_counts))
  
  if (length(missing_markers) > 0) {
    cat("Warning: Missing marker genes:", paste(missing_markers, collapse = ", "), "\n")
  }
  
  if (length(present_markers) == 0) {
    cat("Error: No marker genes found in dataset!\n")
    return(NULL)
  }
  
  # Extract marker gene expression
  marker_expression <- cpm_counts[present_markers, , drop = FALSE]
  
  # Create expression data frame
  marker_data <- data.frame(
    sample_id = rep(colnames(marker_expression), each = length(present_markers)),
    gene = rep(present_markers, ncol(marker_expression)),
    expression = as.vector(marker_expression),
    stringsAsFactors = FALSE
  )
  
  # Add sample metadata
  marker_data <- marker_data %>%
    left_join(sample_metadata, by = "sample_id")
  
  return(marker_data)
}

# Analyze marker genes
marker_data <- analyze_marker_genes(combined_counts, sample_metadata, marker_genes)

if (!is.null(marker_data)) {
  # Create boxplot for marker gene expression
  cat("Creating marker gene expression plots...\n")
  
  marker_boxplot <- ggplot(marker_data, aes(x = cancer_type, y = expression, fill = cancer_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = c("LUAD" = "#E31A1C", "LUSC" = "#1F78B4")) +
    labs(
      title = "Marker Gene Expression by Cancer Type",
      x = "Cancer Type",
      y = "log2(CPM+1)",
      fill = "Cancer Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      strip.text = element_text(size = 12)
    )
  
  ggsave(file.path(plots_dir, "marker_genes_boxplot.png"), marker_boxplot, 
         width = 12, height = 8, dpi = 300)
  
  # Create violin plot
  marker_violin <- ggplot(marker_data, aes(x = cancer_type, y = expression, fill = cancer_type)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = c("LUAD" = "#E31A1C", "LUSC" = "#1F78B4")) +
    labs(
      title = "Marker Gene Expression Distribution by Cancer Type",
      x = "Cancer Type",
      y = "log2(CPM+1)",
      fill = "Cancer Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      strip.text = element_text(size = 12)
    )
  
  ggsave(file.path(plots_dir, "marker_genes_violin.png"), marker_violin, 
         width = 12, height = 8, dpi = 300)
  
  # Statistical analysis of marker genes
  cat("Performing statistical analysis...\n")
  
  # Calculate summary statistics
  marker_summary <- marker_data %>%
    group_by(gene, cancer_type) %>%
    summarise(
      mean_expression = mean(expression),
      median_expression = median(expression),
      sd_expression = sd(expression),
      n_samples = n(),
      .groups = "drop"
    )
  
  # Perform t-tests for each marker gene
  marker_tests <- marker_data %>%
    group_by(gene) %>%
    do(
      ttest_result = t.test(expression ~ cancer_type, data = .)
    ) %>%
    mutate(
      p_value = sapply(ttest_result, function(x) x$p.value),
      t_statistic = sapply(ttest_result, function(x) x$statistic),
      luad_mean = sapply(ttest_result, function(x) x$estimate[1]),
      lusc_mean = sapply(ttest_result, function(x) x$estimate[2])
    ) %>%
    select(-ttest_result)
  
  # Adjust p-values for multiple testing
  marker_tests$p_adjusted <- p.adjust(marker_tests$p_value, method = "BH")
  
  # Save results
  write.csv(marker_data, file = file.path(tables_dir, "marker_gene_expression.csv"), row.names = FALSE)
  write.csv(marker_summary, file = file.path(tables_dir, "marker_gene_summary.csv"), row.names = FALSE)
  write.csv(marker_tests, file = file.path(tables_dir, "marker_gene_tests.csv"), row.names = FALSE)
  
  # Create PCA plot colored by PC1 scores
  cat("Creating PCA plots with PC1 coloring...\n")
  
  # Add PC1 scores to sample metadata
  pc1_data <- cpm_pca$pc_scores %>%
    select(sample_id, PC1, cancer_type)
  
  # Create PCA plot colored by PC1
  pc1_plot <- ggplot(pc1_data, aes(x = PC1, y = 0, color = PC1)) +
    geom_point(size = 3, alpha = 0.7) +
    facet_wrap(~cancer_type) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = "PC1 Scores by Cancer Type",
      x = "PC1 Score",
      y = "",
      color = "PC1 Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    )
  
  ggsave(file.path(plots_dir, "pc1_scores_by_cancer_type.png"), pc1_plot, 
         width = 10, height = 6, dpi = 300)
  
  cat("Marker gene analysis completed successfully!\n")
  cat("Generated plots:\n")
  cat("- marker_genes_boxplot.png\n")
  cat("- marker_genes_violin.png\n")
  cat("- pc1_scores_by_cancer_type.png\n")
  cat("Generated tables:\n")
  cat("- marker_gene_expression.csv\n")
  cat("- marker_gene_summary.csv\n")
  cat("- marker_gene_tests.csv\n")
  
} else {
  cat("Marker gene analysis failed - no marker genes found!\n")
}