#!/usr/bin/env Rscript

# Diagnostic script to examine patient ID mapping issues
suppressMessages({
  library(dplyr)
})

# Load the existing data to inspect actual formats
multiomics_dir <- file.path("data", "multiomics")

# Function to safely load data
safe_load <- function(file_path) {
  if (file.exists(file_path)) {
    load(file_path)
    return(TRUE)
  }
  return(FALSE)
}

# Diagnostic results storage
diagnostics <- list()

# 1. Examine mRNA data
if (safe_load(file.path(multiomics_dir, "mrna_data.RData"))) {
  cat("=== mRNA DATA INSPECTION ===\n")
  cat("Column names (first 10):", paste(head(colnames(mrna_counts), 10), collapse=", "), "\n")
  cat("Sample metadata columns:", paste(colnames(mrna_samples), collapse=", "), "\n")
  
  # Show actual sample names
  cat("First 5 actual sample names:", paste(head(colnames(mrna_counts), 5), collapse=", "), "\n")
  
  # Extract patient IDs using current logic
  mrna_patients <- substr(colnames(mrna_counts), 1, 12)
  cat("Extracted patient IDs (first 10):", paste(head(unique(mrna_patients), 10), collapse=", "), "\n")
  cat("Unique mRNA patients:", length(unique(mrna_patients)), "\n\n")
  
  diagnostics$mrna <- list(
    samples = colnames(mrna_counts),
    patients = mrna_patients,
    unique_patients = unique(mrna_patients)
  )
}

# 2. Examine miRNA data  
if (safe_load(file.path(multiomics_dir, "mirna_data.RData"))) {
  cat("=== miRNA DATA INSPECTION ===\n")
  cat("Column names (first 10):", paste(head(colnames(mirna_counts), 10), collapse=", "), "\n")
  cat("Sample metadata columns:", paste(colnames(mirna_samples), collapse=", "), "\n")
  
  # Show actual sample names
  cat("First 5 actual sample names:", paste(head(colnames(mirna_counts), 5), collapse=", "), "\n")
  
  # Extract patient IDs using current logic
  mirna_patients <- substr(colnames(mirna_counts), 1, 12)
  cat("Extracted patient IDs (first 10):", paste(head(unique(mirna_patients), 10), collapse=", "), "\n")
  cat("Unique miRNA patients:", length(unique(mirna_patients)), "\n\n")
  
  diagnostics$mirna <- list(
    samples = colnames(mirna_counts),
    patients = mirna_patients,
    unique_patients = unique(mirna_patients)
  )
}

# 3. Examine protein data
if (safe_load(file.path(multiomics_dir, "protein_data.RData"))) {
  cat("=== PROTEIN DATA INSPECTION ===\n")
  cat("Column names (first 10):", paste(head(colnames(protein_matrix), 10), collapse=", "), "\n")
  cat("Sample metadata columns:", paste(colnames(protein_samples), collapse=", "), "\n")
  
  # Show actual sample names
  cat("First 5 actual sample names:", paste(head(colnames(protein_matrix), 5), collapse=", "), "\n")
  
  # Extract patient IDs using current logic
  protein_patients <- substr(colnames(protein_matrix), 1, 12)
  cat("Extracted patient IDs (first 10):", paste(head(unique(protein_patients), 10), collapse=", "), "\n")
  cat("Unique protein patients:", length(unique(protein_patients)), "\n\n")
  
  diagnostics$protein <- list(
    samples = colnames(protein_matrix),
    patients = protein_patients,
    unique_patients = unique(protein_patients)
  )
}

# 4. Check subtype data
if (safe_load(file.path(multiomics_dir, "subtype_info.RData"))) {
  cat("=== SUBTYPE DATA INSPECTION ===\n")
  cat("Subtype data columns:", paste(colnames(subtype_info), collapse=", "), "\n")
  cat("First 5 patient IDs:", paste(head(subtype_info$patient_id, 5), collapse=", "), "\n")
  cat("Unique subtypes:", paste(unique(subtype_info$subtype), collapse=", "), "\n\n")
  
  diagnostics$subtype <- list(
    patient_ids = subtype_info$patient_id,
    unique_patients = unique(subtype_info$patient_id),
    subtypes = unique(subtype_info$subtype)
  )
}

# 5. Analyze common patients
if (all(c("mrna", "mirna", "protein", "subtype") %in% names(diagnostics))) {
  cat("=== COMMON PATIENT ANALYSIS ===\n")
  
  # Find common patients
  common_mrna_mirna <- intersect(diagnostics$mrna$unique_patients, diagnostics$mirna$unique_patients)
  common_all <- intersect(intersect(diagnostics$mrna$unique_patients, diagnostics$mirna$unique_patients), 
                         diagnostics$protein$unique_patients)
  common_with_subtype <- intersect(common_all, diagnostics$subtype$unique_patients)
  
  cat("Common between mRNA and miRNA:", length(common_mrna_mirna), "\n")
  cat("Common across all three data types:", length(common_all), "\n")
  cat("Common across all data types with subtype info:", length(common_with_subtype), "\n")
  
  # Show examples of mismatched formats
  cat("\nSample format analysis:\n")
  cat("mRNA format example:", head(diagnostics$mrna$unique_patients, 3), "\n")
  cat("miRNA format example:", head(diagnostics$mirna$unique_patients, 3), "\n")
  cat("protein format example:", head(diagnostics$protein$unique_patients, 3), "\n")
  cat("subtype format example:", head(diagnostics$subtype$unique_patients, 3), "\n")
  
  # Detailed format comparison
  cat("\nFormat patterns:\n")
  cat("mRNA TCGA pattern:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", diagnostics$mrna$unique_patients)), "/", length(diagnostics$mrna$unique_patients), "\n")
  cat("miRNA TCGA pattern:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", diagnostics$mirna$unique_patients)), "/", length(diagnostics$mirna$unique_patients), "\n")
  cat("protein TCGA pattern:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", diagnostics$protein$unique_patients)), "/", length(diagnostics$protein$unique_patients), "\n")
  cat("subtype TCGA pattern:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", diagnostics$subtype$unique_patients)), "/", length(diagnostics$subtype$unique_patients), "\n")
}

# 6. Save diagnostic results
save(diagnostics, file = file.path("logs", "patient_id_diagnostics.RData"))
cat("\nDiagnostic results saved to logs/patient_id_diagnostics.RData\n")