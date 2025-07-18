#!/usr/bin/env Rscript

# Comprehensive Data Inspection Report
# Purpose: Document exact sample ID formats and mapping issues

suppressMessages({
  library(dplyr)
})

# Create detailed inspection report
report_file <- "logs/data_inspection_report.txt"
sink(report_file)

cat(strrep("=", 80), "\n")
cat("TCGA PATIENT ID MAPPING INSPECTION REPORT\n")
cat(strrep("=", 80), "\n\n")

data_dir <- "data/multiomics"

# =============================================================================
# 1. mRNA DATA INSPECTION
# =============================================================================

tryCatch({
  load(file.path(data_dir, "mrna_data.RData"))
  
  cat("1. MRNA DATA:\n")
  cat(strrep("-", 60), "\n")
  cat("Matrix dimensions:", nrow(mrna_counts), "x", ncol(mrna_counts), "\n")
  cat("Column names (first 10):", paste(head(colnames(mrna_counts), 10), collapse=", "), "\n")
  cat("\nData frame head (first 5 samples, first 3 genes):\n")
  print(mrna_counts[1:3, 1:5])
  
  cat("\nPATIENT ID FORMAT ANALYSIS:\n")
  all_patient_ids <- unique(substr(colnames(mrna_counts), 1, 12))
  cat("Total unique patients:", length(all_patient_ids), "\n")
  
  tcga_format <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", all_patient_ids)
  cat("TCGA format (TCGA-XX-XXXX):", sum(tcga_format), "\n")
  cat("First 10 patient IDs:", paste(head(all_patient_ids, 10), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("ERROR loading mRNA data:", e$message, "\n")
})

cat("\n" %>% strrep("=", 80) %>% "\n\n")

# =============================================================================
# 2. miRNA DATA INSPECTION
# =============================================================================

tryCatch({
  load(file.path(data_dir, "mirna_data.RData"))
  
  cat("2. miRNA DATA:\n")
  cat(strrep("-", 60), "\n")
  cat("Matrix dimensions:", nrow(mirna_counts), "x", ncol(mirna_counts), "\n")
  cat("Column names (first 10):", paste(head(colnames(mirna_counts), 10), collapse=", "), "\n")
  cat("\nData frame head (first 5 samples, first 3 miRNAs):\n")
  print(mirna_counts[1:3, 1:5])
  
  cat("\nPATIENT ID FORMAT ANALYSIS:\n")
  all_patient_ids <- unique(substr(colnames(mirna_counts), 1, 12))
  cat("Total unique patients:", length(all_patient_ids), "\n")
  
  tcga_format <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", all_patient_ids)
  uuid_format <- grepl("^[0-9a-f]{8}-", all_patient_ids)
  
  cat("TCGA format (TCGA-XX-XXXX):", sum(tcga_format), "\n")
  cat("UUID format (UUID prefix):", sum(uuid_format), "\n")
  cat("First 10 patient IDs:", paste(head(all_patient_ids, 10), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("ERROR loading miRNA data:", e$message, "\n")
})

cat("\n" %>% strrep("=", 80) %>% "\n\n")

# =============================================================================
# 3. PROTEIN DATA INSPECTION
# =============================================================================

tryCatch({
  load(file.path(data_dir, "protein_data.RData"))
  
  cat("3. PROTEIN DATA:\n")
  cat(strrep("-", 60), "\n")
  cat("Matrix dimensions:", nrow(protein_matrix), "x", ncol(protein_matrix), "\n")
  cat("Column names (first 10):", paste(head(colnames(protein_matrix), 10), collapse=", "), "\n")
  cat("\nData frame head (first 5 samples, first 3 proteins):\n")
  print(protein_matrix[1:3, 1:5])
  
  cat("\nPATIENT ID FORMAT ANALYSIS:\n")
  all_patient_ids <- unique(substr(colnames(protein_matrix), 1, 12))
  cat("Total unique patients:", length(all_patient_ids), "\n")
  
  tcga_format <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", all_patient_ids)
  cat("TCGA format (TCGA-XX-XXXX):", sum(tcga_format), "\n")
  cat("First 10 patient IDs:", paste(head(all_patient_ids, 10), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("ERROR loading protein data:", e$message, "\n")
})

cat("\n" %>% strrep("=", 80) %>% "\n\n")

# =============================================================================
# 4. SUBTYPE DATA INSPECTION
# =============================================================================

tryCatch({
  load(file.path(data_dir, "subtype_info.RData"))
  
  cat("4. SUBTYPE DATA:\n")
  cat(strrep("-", 60), "\n")
  cat("Data frame dimensions:", nrow(subtype_info), "x", ncol(subtype_info), "\n")
  cat("Columns:", paste(colnames(subtype_info), collapse = ", "), "\n")
  cat("\nData head:\n")
  print(head(subtype_info, 5))
  
  cat("\nPATIENT ID FORMAT ANALYSIS:\n")
  all_patient_ids <- unique(subtype_info$patient_id)
  cat("Total unique patients:", length(all_patient_ids), "\n")
  
  tcga_format <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", all_patient_ids)
  cat("TCGA format (TCGA-XX-XXXX):", sum(tcga_format), "\n")
  cat("First 10 patient IDs:", paste(head(all_patient_ids, 10), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("ERROR loading subtype data:", e$message, "\n")
})

cat("\n" %>% strrep("=", 80) %>% "\n\n")

# =============================================================================
# 5. PATIENT OVERLAP ANALYSIS
# =============================================================================

tryCatch({
  # Reload all data for overlap analysis
  load(file.path(data_dir, "mrna_data.RData"))
  load(file.path(data_dir, "mirna_data.RData"))
  load(file.path(data_dir, "protein_data.RData"))
  load(file.path(data_dir, "subtype_info.RData"))
  
  cat("5. PATIENT OVERLAP ANALYSIS:\n")
  cat(strrep("-", 60), "\n")
  
  # Extract patient IDs consistently
  mrna_patients <- unique(substr(colnames(mrna_counts), 1, 12))
  mirna_patients <- unique(substr(colnames(mirna_counts), 1, 12))
  protein_patients <- unique(substr(colnames(protein_matrix), 1, 12))
  subtype_patients <- unique(subtype_info$patient_id)
  
  cat("Dataset Sizes:\n")
  cat("  mRNA:  ", length(mrna_patients), "patients\n")
  cat("  miRNA: ", length(mirna_patients), "patients\n")
  cat("  protein:", length(protein_patients), "patients\n")
  cat("  subtype:", length(subtype_patients), "patients\n\n")
  
  # Calculate overlaps
  mrna_mirna_overlap <- length(intersect(mrna_patients, mirna_patients))
  mrna_protein_overlap <- length(intersect(mrna_patients, protein_patients))
  mirna_protein_overlap <- length(intersect(mirna_patients, protein_patients))
  
  all_three_overlap <- length(Reduce(intersect, list(mrna_patients, mirna_patients, protein_patients)))
  all_with_subtype <- length(Reduce(intersect, list(mrna_patients, mirna_patients, protein_patients, subtype_patients)))
  
  cat("Overlap Analysis:\n")
  cat("  mRNA ∩ miRNA:", mrna_mirna_overlap, "patients\n")
  cat("  mRNA ∩ protein:", mrna_protein_overlap, "patients\n")
  cat("  miRNA ∩ protein:", mirna_protein_overlap, "patients\n")
  cat("  All three datasets:", all_three_overlap, "patients\n")
  cat("  All datasets + subtype:", all_with_subtype, "patients\n\n")
  
  if (all_three_overlap == 0) {
    cat("🔴 CRITICAL: NO OVERLAP FOUND BETWEEN DATASETS!\n\n")
    
    # Show format comparison
    cat("Format Comparison:\n")
    cat("  mRNA IDs:", head(mrna_patients, 3), "...\n")
    cat("  miRNA IDs:", head(mirna_patients, 3), "...\n")
    cat("  protein IDs:", head(protein_patients, 3), "...\n")
    cat("  subtype IDs:", head(subtype_patients, 3), "...\n\n")
    
    # Detailed format analysis
    formats <- data.frame(
      Dataset = c("mRNA", "miRNA", "protein", "subtype"),
      TCGA_format_count = c(
        sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", mrna_patients)),
        sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", mirna_patients)),
        sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", protein_patients)),
        sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", subtype_patients))
      ),
      Total = c(length(mrna_patients), length(mirna_patients), 
                length(protein_patients), length(subtype_patients))
    )
    
    formats$TCGA_percentage <- sprintf("%.1f%%", formats$TCGA_format_count / formats$Total * 100)
    
    cat("Detailed Format Analysis:\n")
    print(formats, row.names = FALSE)
  } else {
    cat("✅ Common patients found:", all_three_overlap, "\n")
  }
  
}, error = function(e) {
  cat("ERROR in overlap analysis:", e$message, "\n")
})

cat("\n" %>% strrep("=", 80) %>% "\n")
cat("END OF REPORT\n")
cat("=" %>% strrep("=", 80) %>% "\n")

sink()

cat("Data inspection report saved to:", report_file, "\n")

# Also create a quick summary
summary_file <- "logs/patient_mapping_summary.txt"
sink(summary_file)

cat("PATIENT MAPPING SUMMARY\n")
cat("======================\n\n")

tryCatch({
  load(file.path(data_dir, "mrna_data.RData"))
  load(file.path(data_dir, "mirna_data.RData"))
  load(file.path(data_dir, "protein_data.RData"))
  load(file.path(data_dir, "subtype_info.RData"))
  
  mrna_patients <- unique(substr(colnames(mrna_counts), 1, 12))
  mirna_patients <- unique(substr(colnames(mirna_counts), 1, 12))
  protein_patients <- unique(substr(colnames(protein_matrix), 1, 12))
  subtype_patients <- unique(subtype_info$patient_id)
  
  cat("Dataset Sizes:\n")
  cat("  mRNA:  ", length(mrna_patients), "patients\n")
  cat("  miRNA: ", length(mirna_patients), "patients\n")
  cat("  protein:", length(protein_patients), "patients\n")
  cat("  subtype:", length(subtype_patients), "patients\n\n")
  
  common_all <- length(Reduce(intersect, list(mrna_patients, mirna_patients, protein_patients)))
  
  if (common_all == 0) {
    cat("❌ ISSUE CONFIRMED: Zero common patients across datasets\n\n")
    cat("Root Cause: miRNA uses UUID format instead of TCGA format\n")
    cat("Solution: Fix miRNA sample ID mapping\n")
  } else {
    cat("✅ Common patients found:", common_all, "\n")
  }
  
}, error = function(e) {
  cat("Error in summary:", e$message, "\n")
})

sink()

cat("Summary saved to:", summary_file, "\n")