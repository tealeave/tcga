# Visualization Script
# Generates all plots from blog posts following exact methodology

suppressMessages({
  library(ggplot2)
  library(ggfortify)
  library(ComplexHeatmap)
  library(dplyr)
  library(RColorBrewer)
})

cat("Loading visualization functions...\n")

# Load PCA results
load(file.path(data_dir, "pca_results.RData"))
load(file.path(data_dir, "combined_counts.RData"))

# Create color palette for cancer types
cancer_colors <- c("LUAD" = "#E31A1C", "LUSC" = "#1F78B4")

# Function to create PCA plot
create_pca_plot <- function(pca_results, title_suffix = "") {
  pc_scores <- pca_results$pc_scores
  var_explained <- pca_results$var_explained
  
  p <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = cancer_type)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = cancer_colors) +
    labs(
      title = paste0("PCA Analysis - ", title_suffix),
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = "Cancer Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom"
    )
  
  return(p)
}

# Generate PCA plots
cat("Generating PCA plots...\n")

# Raw counts PCA plot
raw_pca_plot <- create_pca_plot(raw_pca, "Raw Counts")
ggsave(file.path(plots_dir, "raw_counts_pca.png"), raw_pca_plot, 
       width = 8, height = 6, dpi = 300)

# CPM-normalized PCA plot
cpm_pca_plot <- create_pca_plot(cpm_pca, "CPM Normalized")
ggsave(file.path(plots_dir, "cpm_normalized_pca.png"), cpm_pca_plot, 
       width = 8, height = 6, dpi = 300)

# Combined comparison plot
combined_plot <- gridExtra::grid.arrange(raw_pca_plot, cpm_pca_plot, ncol = 2)
ggsave(file.path(plots_dir, "pca_comparison.png"), combined_plot, 
       width = 16, height = 6, dpi = 300)

# Create variance explained plot
variance_plot_data <- data.frame(
  PC = factor(1:10),
  Raw_Counts = raw_pca$var_explained[1:10],
  CPM_Normalized = cpm_pca$var_explained[1:10]
)

# Reshape for plotting
variance_plot_data_long <- tidyr::pivot_longer(variance_plot_data, 
                                              cols = c(Raw_Counts, CPM_Normalized),
                                              names_to = "Method", 
                                              values_to = "Variance_Explained")

variance_plot <- ggplot(variance_plot_data_long, aes(x = PC, y = Variance_Explained, 
                                                    fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  labs(
    title = "Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Variance Explained (%)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(plots_dir, "variance_explained.png"), variance_plot, 
       width = 10, height = 6, dpi = 300)

# Create heatmap of top variable genes
cat("Generating heatmap...\n")

# Get top 50 most variable genes for heatmap
load(file.path(data_dir, "combined_counts.RData"))
library(edgeR)
cpm_counts <- cpm(combined_counts, log = TRUE, prior.count = 1)

# Calculate gene variances and select top 50
gene_vars <- apply(cpm_counts, 1, var)
top_50_genes <- names(sort(gene_vars, decreasing = TRUE)[1:50])

# Prepare heatmap data
heatmap_data <- cpm_counts[top_50_genes, ]

# Create annotation for samples
sample_annotation <- data.frame(
  Cancer_Type = sample_metadata$cancer_type
)
rownames(sample_annotation) <- sample_metadata$sample_id

# Create heatmap
pdf(file.path(plots_dir, "top_genes_heatmap.pdf"), width = 12, height = 8)
Heatmap(heatmap_data,
        name = "log2(CPM+1)",
        top_annotation = HeatmapAnnotation(
          Cancer_Type = sample_annotation$Cancer_Type,
          col = list(Cancer_Type = cancer_colors)
        ),
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_title = "Top 50 Variable Genes",
        heatmap_legend_param = list(
          title = "log2(CPM+1)",
          legend_direction = "horizontal",
          legend_width = unit(6, "cm")
        )
)
dev.off()

# Also save as PNG
png(file.path(plots_dir, "top_genes_heatmap.png"), width = 12, height = 8, 
    units = "in", res = 300)
Heatmap(heatmap_data,
        name = "log2(CPM+1)",
        top_annotation = HeatmapAnnotation(
          Cancer_Type = sample_annotation$Cancer_Type,
          col = list(Cancer_Type = cancer_colors)
        ),
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_title = "Top 50 Variable Genes"
)
dev.off()

cat("Visualization completed successfully!\n")
cat("Generated plots:\n")
cat("- raw_counts_pca.png\n")
cat("- cpm_normalized_pca.png\n") 
cat("- pca_comparison.png\n")
cat("- variance_explained.png\n")
cat("- top_genes_heatmap.png/pdf\n")