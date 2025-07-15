#!/usr/bin/env Rscript

# TCGA RNA-seq PCA Analysis Pipeline
# Main driver script for TCGA LUAD/LUSC analysis following blog methodology

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
cat("TCGA RNA-seq PCA Analysis Pipeline\n")
cat("========================================\n\n")

# Step 1: Data Download
cat("Step 1: Downloading TCGA data...\n")
source(file.path(src_dir, "01_data_download.R"))

# Step 2: Data Preprocessing
cat("Step 2: Preprocessing data...\n")
source(file.path(src_dir, "02_data_preprocessing.R"))

# Step 3: PCA Analysis
cat("Step 3: Running PCA analysis...\n")
source(file.path(src_dir, "03_pca_analysis.R"))

# Step 4: Visualization
cat("Step 4: Generating visualizations...\n")
source(file.path(src_dir, "04_visualization.R"))

# Step 5: Marker Analysis
cat("Step 5: Running marker gene analysis...\n")
source(file.path(src_dir, "05_marker_analysis.R"))

cat("\n========================================\n")
cat("Pipeline completed successfully!\n")
cat("Results saved to:\n")
cat("- Plots: ", plots_dir, "\n")
cat("- Tables: ", tables_dir, "\n")
cat("========================================\n")