#!/usr/bin/env Rscript

# TCGA Multi-omics DIABLO Analysis Pipeline
# Main driver script for TCGA breast cancer multi-omics analysis using DIABLO

# Load required libraries
suppressMessages({
  library(renv)
  renv::activate()
})

# Set up paths
project_root <- getwd()
src_dir <- file.path(project_root, "src")
data_dir <- file.path(project_root, "data")
results_dir <- file.path(project_root, "results")
plots_dir <- file.path(results_dir, "plots")
tables_dir <- file.path(results_dir, "tables")

# Create output directories if they don't exist
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Pipeline execution
cat("========================================\n")
cat("TCGA Multi-omics DIABLO Analysis Pipeline\n")
cat("========================================\n\n")

# Record start time
start_time <- Sys.time()

# Step 1: Multi-omics Data Download
cat("Step 1: Downloading multi-omics TCGA data...\n")
cat("  - Downloading breast cancer mRNA, miRNA, and protein data\n")
cat("  - Creating training and test sets\n")
cat("  - Extracting breast cancer subtype information\n")
source(file.path(src_dir, "06_multiomics_download.R"))

# Step 2: Multi-omics Data Preprocessing
cat("\nStep 2: Preprocessing multi-omics data...\n")
cat("  - Standardizing data matrices across omics types\n")
cat("  - Handling sample matching across data types\n")
cat("  - Creating design matrix for DIABLO\n")
source(file.path(src_dir, "07_multiomics_preprocessing.R"))

# Step 3: Block sPLS-DA Analysis
cat("\nStep 3: Running Block sPLS-DA analysis...\n")
cat("  - Tuning number of components\n")
cat("  - Optimizing variable selection\n")
cat("  - Building final multiblock model\n")
source(file.path(src_dir, "08_block_splsda.R"))

# Step 4: DIABLO Analysis
cat("\nStep 4: Running DIABLO analysis...\n")
cat("  - Implementing full DIABLO workflow\n")
cat("  - Handling missing data blocks\n")
cat("  - Generating weighted and majority vote predictions\n")
source(file.path(src_dir, "09_diablo_analysis.R"))

# Step 5: Multi-omics Visualization
cat("\nStep 5: Generating multi-omics visualizations...\n")
cat("  - Creating sample plots and network visualizations\n")
cat("  - Generating clustered image maps\n")
cat("  - Producing performance and biomarker plots\n")
source(file.path(src_dir, "10_multiomics_visualization.R"))

# Calculate runtime
end_time <- Sys.time()
runtime <- end_time - start_time

# Pipeline summary
cat("\n========================================\n")
cat("Pipeline completed successfully!\n")
cat("========================================\n")
cat("Runtime: ", format(runtime), "\n")
cat("Results saved to:\n")
cat("- Multi-omics data: ", file.path(data_dir, "multiomics"), "\n")
cat("- Plots: ", plots_dir, "\n")
cat("- Tables: ", tables_dir, "\n")

# Load final results for summary
multiomics_dir <- file.path(data_dir, "multiomics")
if (file.exists(file.path(multiomics_dir, "diablo_results.RData"))) {
  load(file.path(multiomics_dir, "diablo_results.RData"))
  
  cat("\n=== Analysis Summary ===\n")
  
  # Data summary
  cat("Data Summary:\n")
  load(file.path(multiomics_dir, "preprocessing_summary.RData"))
  
  cat("  Training set:\n")
  for (block in names(train_summary)) {
    if (block != "outcome" && block != "dataset") {
      cat("    ", block, ": ", train_summary[[block]]$n_features, " features x ", 
          train_summary[[block]]$n_samples, " samples\n")
    }
  }
  cat("    Outcome distribution: ", paste(names(train_summary$outcome), 
                                       train_summary$outcome, sep = "=", collapse = ", "), "\n")
  
  cat("  Test set:\n")
  for (block in names(test_summary)) {
    if (block != "outcome" && block != "dataset") {
      cat("    ", block, ": ", test_summary[[block]]$n_features, " features x ", 
          test_summary[[block]]$n_samples, " samples\n")
    }
  }
  cat("    Outcome distribution: ", paste(names(test_summary$outcome), 
                                       test_summary$outcome, sep = "=", collapse = ", "), "\n")
  
  # Model summary
  cat("\nDIABLO Model Summary:\n")
  cat("  Optimal components: ", diablo_final_results$main_results$optimal_ncomp, "\n")
  cat("  keepX parameters:\n")
  for (block_name in names(diablo_final_results$main_results$optimal_keepx)) {
    cat("    ", block_name, ": ", paste(diablo_final_results$main_results$optimal_keepx[[block_name]], collapse = ", "), "\n")
  }
  
  # Performance summary
  cat("\nPerformance Summary:\n")
  cat("  Cross-validation BER: ", round(diablo_final_results$main_results$performance$error.rate$BER[diablo_final_results$main_results$optimal_ncomp], 3), "\n")
  
  if (!is.null(diablo_final_results$predictions) && !is.null(diablo_final_results$predictions$evaluation)) {
    cat("  Test set predictions:\n")
    cat("    Centroids distance: ", round(diablo_final_results$predictions$evaluation$accuracy_centroids * 100, 2), "%\n")
    cat("    Majority vote: ", round(diablo_final_results$predictions$evaluation$accuracy_majority * 100, 2), "%\n")
    cat("    Weighted vote: ", round(diablo_final_results$predictions$evaluation$accuracy_weighted * 100, 2), "%\n")
  }
  
  # Top biomarkers
  cat("\nTop Biomarkers (per data block):\n")
  for (block_name in names(diablo_final_results$biomarkers)) {
    top_biomarkers <- diablo_final_results$biomarkers[[block_name]]$biomarkers[1:min(5, length(diablo_final_results$biomarkers[[block_name]]$biomarkers))]
    cat("  ", block_name, ": ", paste(top_biomarkers, collapse = ", "), "\n")
  }
  
  # Generated outputs
  cat("\nGenerated Outputs:\n")
  
  # Count plots
  plot_files <- list.files(plots_dir, pattern = "*.png", full.names = FALSE)
  diablo_plots <- sum(grepl("diablo", plot_files))
  block_plots <- sum(grepl("block", plot_files))
  
  cat("  Plots generated:\n")
  cat("    DIABLO plots: ", diablo_plots, "\n")
  cat("    Block sPLS-DA plots: ", block_plots, "\n")
  cat("    Total plots: ", length(plot_files), "\n")
  
  # Count tables
  table_files <- list.files(tables_dir, pattern = "*.csv", full.names = FALSE)
  diablo_tables <- sum(grepl("diablo", table_files))
  block_tables <- sum(grepl("block", table_files))
  
  cat("  Tables generated:\n")
  cat("    DIABLO tables: ", diablo_tables, "\n")
  cat("    Block sPLS-DA tables: ", block_tables, "\n")
  cat("    Total tables: ", length(table_files), "\n")
  
} else {
  cat("Warning: DIABLO results not found. Pipeline may have failed.\n")
}

cat("\n========================================\n")
cat("For detailed analysis, examine the files in:\n")
cat("- ", file.path(data_dir, "multiomics"), " (analysis results)\n")
cat("- ", plots_dir, " (visualizations)\n")
cat("- ", tables_dir, " (summary tables)\n")
cat("========================================\n")

# Optional: Create a quick report
create_quick_report <- function() {
  report_content <- paste0(
    "# TCGA Multi-omics DIABLO Analysis Report\n\n",
    "**Analysis completed on:** ", Sys.Date(), "\n",
    "**Runtime:** ", format(runtime), "\n\n",
    "## Pipeline Overview\n",
    "This analysis implemented the DIABLO (Data Integration Analysis for Biomarker discovery using Latent variable approaches) methodology on TCGA breast cancer data.\n\n",
    "### Data Types Analyzed\n",
    "- mRNA expression data\n",
    "- miRNA expression data\n",
    "- Protein expression data\n\n",
    "### Breast Cancer Subtypes\n",
    "- Basal\n",
    "- Her2\n",
    "- LumA\n\n",
    "### Key Results\n",
    "- Model successfully trained on multi-omics data\n",
    "- Biomarker signatures identified for each data type\n",
    "- Cross-validation and test set performance evaluated\n",
    "- Comprehensive visualizations generated\n\n",
    "### Files Generated\n",
    "- Analysis results: `data/multiomics/`\n",
    "- Visualizations: `results/plots/`\n",
    "- Summary tables: `results/tables/`\n\n",
    "For detailed results, examine the generated files and visualizations.\n"
  )
  
  writeLines(report_content, file.path(results_dir, "multiomics_analysis_report.md"))
  cat("Quick report saved to:", file.path(results_dir, "multiomics_analysis_report.md"), "\n")
}

# Create quick report
create_quick_report()

cat("\nMulti-omics DIABLO analysis pipeline completed successfully!\n")