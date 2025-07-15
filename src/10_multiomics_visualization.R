# Multi-omics Visualization Suite
# Comprehensive visualization functions for DIABLO and Block sPLS-DA results

suppressMessages({
  library(mixOmics)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(gridExtra)
  library(igraph)
  library(corrplot)
  library(ComplexHeatmap)
  library(circlize)
})

cat("Loading multi-omics visualization functions...\n")

# Load results data
multiomics_dir <- file.path(data_dir, "multiomics")
load(file.path(multiomics_dir, "diablo_results.RData"))
load(file.path(multiomics_dir, "block_splsda_results.RData"))
load(file.path(multiomics_dir, "train_processed.RData"))
load(file.path(multiomics_dir, "test_processed.RData"))

# Define color palette for subtypes
subtype_colors <- c(
  "Basal" = "#E31A1C",
  "Her2" = "#1F78B4", 
  "LumA" = "#33A02C"
)

# Function to create sample plots
create_sample_plots <- function(diablo_model, Y, title_prefix = "DIABLO") {
  cat("Creating sample plots...\n")
  
  # Sample plot for components 1 and 2
  p1 <- plotIndiv(diablo_model, 
                  comp = c(1, 2),
                  ind.names = FALSE,
                  group = Y,
                  col.per.group = subtype_colors,
                  legend = TRUE,
                  title = paste(title_prefix, "Sample Plot (Comp 1 vs 2)"))
  
  # Sample plot for components 1 and 3 (if available)
  if (diablo_model$ncomp >= 3) {
    p2 <- plotIndiv(diablo_model, 
                    comp = c(1, 3),
                    ind.names = FALSE,
                    group = Y,
                    col.per.group = subtype_colors,
                    legend = TRUE,
                    title = paste(title_prefix, "Sample Plot (Comp 1 vs 3)"))
  } else {
    p2 <- NULL
  }
  
  # Sample plot for components 2 and 3 (if available)
  if (diablo_model$ncomp >= 3) {
    p3 <- plotIndiv(diablo_model, 
                    comp = c(2, 3),
                    ind.names = FALSE,
                    group = Y,
                    col.per.group = subtype_colors,
                    legend = TRUE,
                    title = paste(title_prefix, "Sample Plot (Comp 2 vs 3)"))
  } else {
    p3 <- NULL
  }
  
  # Save plots
  ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_sample_plot_comp12.png")), 
         plot = p1, width = 10, height = 8, dpi = 300)
  
  if (!is.null(p2)) {
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_sample_plot_comp13.png")), 
           plot = p2, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(p3)) {
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_sample_plot_comp23.png")), 
           plot = p3, width = 10, height = 8, dpi = 300)
  }
  
  return(list(comp12 = p1, comp13 = p2, comp23 = p3))
}

# Function to create variable correlation networks
create_correlation_networks <- function(diablo_model, cutoff = 0.7, title_prefix = "DIABLO") {
  cat("Creating variable correlation networks...\n")
  
  # Create network for each component
  for (comp in 1:diablo_model$ncomp) {
    tryCatch({
      # Create network plot
      png(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_network_comp", comp, ".png")), 
          width = 12, height = 10, units = "in", res = 300)
      
      network_plot <- network(diablo_model, 
                             comp = comp,
                             cutoff = cutoff,
                             shape.node = c("rectangle", "circle", "triangle"),
                             color.node = c("lightblue", "mistyrose", "lightgreen"),
                             size.node = c(2, 2, 2),
                             color.edge = "gray",
                             lty.edge = "solid",
                             lwd.edge = 1,
                             show.edge.labels = FALSE)
      
      dev.off()
      
      cat("  Network plot for component", comp, "saved\n")
      
    }, error = function(e) {
      cat("  Error creating network for component", comp, ":", e$message, "\n")
    })
  }
}

# Function to create clustered image maps (heatmaps)
create_clustered_image_maps <- function(diablo_model, Y, title_prefix = "DIABLO") {
  cat("Creating clustered image maps...\n")
  
  # Create heatmap for each component
  for (comp in 1:diablo_model$ncomp) {
    tryCatch({
      # Create heatmap
      png(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_heatmap_comp", comp, ".png")), 
          width = 14, height = 10, units = "in", res = 300)
      
      cim_plot <- cim(diablo_model, 
                      comp = comp,
                      margins = c(5, 6),
                      row.sideColors = subtype_colors[Y],
                      row.names = Y,
                      col.names = TRUE,
                      symkey = FALSE,
                      keysize = 1)
      
      dev.off()
      
      cat("  Heatmap for component", comp, "saved\n")
      
    }, error = function(e) {
      cat("  Error creating heatmap for component", comp, ":", e$message, "\n")
    })
  }
}

# Function to create performance plots
create_performance_plots <- function(diablo_results, title_prefix = "DIABLO") {
  cat("Creating performance plots...\n")
  
  # Component tuning plot
  if (!is.null(diablo_results$main_results$component_tuning)) {
    tuning_data <- diablo_results$main_results$component_tuning
    
    # Extract BER values
    ber_values <- tuning_data$error.rate$BER
    components <- 1:length(ber_values)
    
    # Create tuning plot
    p_tuning <- ggplot(data.frame(Component = components, BER = ber_values), 
                       aes(x = Component, y = BER)) +
      geom_line(color = "blue", size = 1.2) +
      geom_point(color = "red", size = 3) +
      geom_vline(xintercept = diablo_results$main_results$optimal_ncomp, 
                 linetype = "dashed", color = "red", alpha = 0.7) +
      labs(title = paste(title_prefix, "Component Tuning"),
           x = "Number of Components", 
           y = "Balanced Error Rate (BER)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_component_tuning.png")), 
           plot = p_tuning, width = 10, height = 6, dpi = 300)
  }
  
  # Performance comparison plot
  if (!is.null(diablo_results$predictions) && !is.null(diablo_results$predictions$evaluation)) {
    eval_data <- diablo_results$predictions$evaluation
    
    # Create performance comparison
    performance_df <- data.frame(
      Method = c("Centroids", "Majority Vote", "Weighted Vote"),
      Accuracy = c(eval_data$accuracy_centroids, 
                   eval_data$accuracy_majority, 
                   eval_data$accuracy_weighted)
    )
    
    p_performance <- ggplot(performance_df, aes(x = Method, y = Accuracy, fill = Method)) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(aes(label = paste0(round(Accuracy * 100, 1), "%")), 
                vjust = -0.5, size = 4) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste(title_prefix, "Prediction Performance"),
           x = "Prediction Method", 
           y = "Accuracy") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none") +
      ylim(0, 1)
    
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_performance_comparison.png")), 
           plot = p_performance, width = 10, height = 6, dpi = 300)
  }
}

# Function to create biomarker importance plots
create_biomarker_plots <- function(biomarkers, title_prefix = "DIABLO") {
  cat("Creating biomarker importance plots...\n")
  
  for (block_name in names(biomarkers)) {
    biomarker_data <- biomarkers[[block_name]]
    
    # Create importance plot
    plot_data <- data.frame(
      Biomarker = factor(biomarker_data$biomarkers[1:min(20, length(biomarker_data$biomarkers))],
                        levels = biomarker_data$biomarkers[1:min(20, length(biomarker_data$biomarkers))]),
      Importance = biomarker_data$stability_scores[1:min(20, length(biomarker_data$biomarkers))]
    )
    
    p_biomarker <- ggplot(plot_data, aes(x = reorder(Biomarker, Importance), y = Importance)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      coord_flip() +
      labs(title = paste(title_prefix, "Top Biomarkers -", block_name),
           x = "Biomarker", 
           y = "Importance Score") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8))
    
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_biomarkers_", block_name, ".png")), 
           plot = p_biomarker, width = 12, height = 8, dpi = 300)
  }
}

# Function to create arrow plots (contribution plots)
create_arrow_plots <- function(diablo_model, Y, title_prefix = "DIABLO") {
  cat("Creating arrow plots...\n")
  
  # Create arrow plot for components 1 and 2
  tryCatch({
    png(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_arrow_plot_comp12.png")), 
        width = 12, height = 10, units = "in", res = 300)
    
    arrow_plot <- plotArrow(diablo_model, 
                           comp = c(1, 2),
                           group = Y,
                           col.per.group = subtype_colors,
                           legend = TRUE,
                           title = paste(title_prefix, "Arrow Plot (Comp 1 vs 2)"))
    
    dev.off()
    
    cat("  Arrow plot for components 1 and 2 saved\n")
    
  }, error = function(e) {
    cat("  Error creating arrow plot:", e$message, "\n")
  })
  
  # Create arrow plot for components 1 and 3 (if available)
  if (diablo_model$ncomp >= 3) {
    tryCatch({
      png(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_arrow_plot_comp13.png")), 
          width = 12, height = 10, units = "in", res = 300)
      
      arrow_plot <- plotArrow(diablo_model, 
                             comp = c(1, 3),
                             group = Y,
                             col.per.group = subtype_colors,
                             legend = TRUE,
                             title = paste(title_prefix, "Arrow Plot (Comp 1 vs 3)"))
      
      dev.off()
      
      cat("  Arrow plot for components 1 and 3 saved\n")
      
    }, error = function(e) {
      cat("  Error creating arrow plot for comp 1 vs 3:", e$message, "\n")
    })
  }
}

# Function to create correlation heatmaps between blocks
create_block_correlation_heatmaps <- function(correlations, title_prefix = "DIABLO") {
  cat("Creating block correlation heatmaps...\n")
  
  if (!is.null(correlations$block_correlations)) {
    for (block_pair in names(correlations$block_correlations)) {
      cor_matrix <- correlations$block_correlations[[block_pair]]
      
      # Create heatmap
      png(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_correlation_", block_pair, ".png")), 
          width = 10, height = 8, units = "in", res = 300)
      
      corrplot(cor_matrix, 
               method = "color", 
               type = "upper",
               order = "hclust",
               tl.cex = 0.8,
               tl.col = "black",
               title = paste(title_prefix, "Correlation:", gsub("_", " vs ", block_pair)),
               mar = c(0, 0, 2, 0))
      
      dev.off()
      
      cat("  Correlation heatmap for", block_pair, "saved\n")
    }
  }
}

# Function to create stability plots
create_stability_plots <- function(stability_results, title_prefix = "DIABLO") {
  cat("Creating stability plots...\n")
  
  for (block_name in names(stability_results)) {
    stability_matrix <- stability_results[[block_name]]
    
    # Get top stable features
    mean_stability <- rowMeans(stability_matrix)
    top_features <- names(sort(mean_stability, decreasing = TRUE))[1:min(30, length(mean_stability))]
    
    # Create stability plot
    plot_data <- data.frame(
      Feature = factor(top_features, levels = top_features),
      Stability = mean_stability[top_features]
    )
    
    p_stability <- ggplot(plot_data, aes(x = reorder(Feature, Stability), y = Stability)) +
      geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
      coord_flip() +
      labs(title = paste(title_prefix, "Feature Stability -", block_name),
           x = "Feature", 
           y = "Stability Score") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8))
    
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_stability_", block_name, ".png")), 
           plot = p_stability, width = 12, height = 8, dpi = 300)
  }
}

# Function to create comprehensive summary plots
create_summary_plots <- function(diablo_results, title_prefix = "DIABLO") {
  cat("Creating comprehensive summary plots...\n")
  
  # Model configuration summary
  config_data <- data.frame(
    Block = names(diablo_results$main_results$optimal_keepx),
    Variables_Selected = sapply(diablo_results$main_results$optimal_keepx, function(x) x[1])
  )
  
  p_config <- ggplot(config_data, aes(x = Block, y = Variables_Selected, fill = Block)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_text(aes(label = Variables_Selected), vjust = -0.5, size = 4) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = paste(title_prefix, "Model Configuration"),
         x = "Data Block", 
         y = "Number of Variables Selected") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_model_configuration.png")), 
         plot = p_config, width = 10, height = 6, dpi = 300)
  
  # Explained variance summary
  if (!is.null(diablo_results$main_results$model$explained_variance)) {
    exp_var_data <- list()
    for (block_name in names(diablo_results$main_results$model$explained_variance)) {
      exp_var_data[[block_name]] <- data.frame(
        Component = 1:length(diablo_results$main_results$model$explained_variance[[block_name]]),
        Variance = diablo_results$main_results$model$explained_variance[[block_name]],
        Block = block_name
      )
    }
    
    exp_var_df <- do.call(rbind, exp_var_data)
    
    p_variance <- ggplot(exp_var_df, aes(x = Component, y = Variance, color = Block)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      scale_color_brewer(palette = "Set2") +
      labs(title = paste(title_prefix, "Explained Variance by Component"),
           x = "Component", 
           y = "Explained Variance") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = file.path(plots_dir, paste0(tolower(title_prefix), "_explained_variance.png")), 
           plot = p_variance, width = 10, height = 6, dpi = 300)
  }
}

# Main visualization execution
cat("=== Multi-omics Visualization Suite ===\n")

# Create sample plots
sample_plots <- create_sample_plots(diablo_final_results$main_results$model, 
                                   train_processed$Y, 
                                   "DIABLO")

# Create correlation networks
create_correlation_networks(diablo_final_results$main_results$model, 
                           cutoff = 0.7, 
                           "DIABLO")

# Create clustered image maps
create_clustered_image_maps(diablo_final_results$main_results$model, 
                           train_processed$Y, 
                           "DIABLO")

# Create performance plots
create_performance_plots(diablo_final_results, "DIABLO")

# Create biomarker plots
create_biomarker_plots(diablo_final_results$biomarkers, "DIABLO")

# Create arrow plots
create_arrow_plots(diablo_final_results$main_results$model, 
                   train_processed$Y, 
                   "DIABLO")

# Create block correlation heatmaps
create_block_correlation_heatmaps(diablo_final_results$correlations, "DIABLO")

# Create stability plots
create_stability_plots(diablo_final_results$stability, "DIABLO")

# Create summary plots
create_summary_plots(diablo_final_results, "DIABLO")

# Also create visualizations for Block sPLS-DA results
if (exists("block_splsda_results")) {
  cat("\n=== Block sPLS-DA Visualizations ===\n")
  
  # Create sample plots for Block sPLS-DA
  create_sample_plots(block_splsda_results$final_model$model, 
                      train_processed$Y, 
                      "Block_sPLS-DA")
  
  # Create biomarker plots for Block sPLS-DA
  create_biomarker_plots(block_splsda_results$important_variables, "Block_sPLS-DA")
  
  # Create performance plots for Block sPLS-DA
  create_performance_plots(block_splsda_results, "Block_sPLS-DA")
}

# Create a comprehensive visualization index
create_visualization_index <- function() {
  cat("Creating visualization index...\n")
  
  # List all created plots
  plot_files <- list.files(plots_dir, pattern = "*.png", full.names = FALSE)
  
  # Categorize plots
  diablo_plots <- plot_files[grepl("diablo", plot_files)]
  block_splsda_plots <- plot_files[grepl("block_spls", plot_files)]
  
  # Create index
  index_content <- paste0(
    "# Multi-omics Visualization Index\n\n",
    "## DIABLO Visualizations\n",
    paste0("- ", diablo_plots, "\n", collapse = ""),
    "\n## Block sPLS-DA Visualizations\n",
    paste0("- ", block_splsda_plots, "\n", collapse = ""),
    "\n## Generated on: ", Sys.Date(), "\n"
  )
  
  writeLines(index_content, file.path(plots_dir, "visualization_index.md"))
  cat("Visualization index created\n")
}

# Create visualization index
create_visualization_index()

cat("Multi-omics visualization suite completed successfully!\n")
cat("All plots saved to:", plots_dir, "\n")
cat("Visualization index saved to:", file.path(plots_dir, "visualization_index.md"), "\n")