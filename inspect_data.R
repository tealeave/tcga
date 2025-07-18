#!/usr/bin/env Rscript

# Quick data inspection script
suppressMessages(library(dplyr))

# Set working directory
data_dir <- "data/multiomics"

# Function to inspect data
cat("=== DATA INSPECTION RESULTS ===\n\n")

# 1. mRNA Data
tryCatch({
    load(file.path(data_dir, "mrna_data.RData"))
    cat("=== mRNA DATA ===\n")
    cat("Sample names (first 5):", paste(head(colnames(mrna_counts), 5), collapse=", "), "\n")
    cat("Patient IDs from substr(1,12):", paste(head(unique(substr(colnames(mrna_counts), 1, 12)), 5), collapse=", "), "\n")
    cat("Total samples:", ncol(mrna_counts), "\n")
    cat("Unique patients:", length(unique(substr(colnames(mrna_counts), 1, 12))), "\n\n")
}, error = function(e) cat("mRNA data error:", e$message, "\n\n"))

# 2. miRNA Data
tryCatch({
    load(file.path(data_dir, "mirna_data.RData"))
    cat("=== miRNA DATA ===\n")
    cat("Sample names (first 5):", paste(head(colnames(mirna_counts), 5), collapse=", "), "\n")
    cat("Patient IDs from substr(1,12):", paste(head(unique(substr(colnames(mirna_counts), 1, 12)), 5), collapse=", "), "\n")
    cat("Total samples:", ncol(mirna_counts), "\n")
    cat("Unique patients:", length(unique(substr(colnames(mirna_counts), 1, 12))), "\n\n")
    
    # Check if samples have TCGA format
    samples <- colnames(mirna_counts)
    tcga_pattern <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", samples)
    cat("TCGA format samples:", sum(tcga_pattern), "/", length(samples), "\n")
    cat("Non-TCGA format samples:", head(samples[!tcga_pattern], 5), "\n\n")
}, error = function(e) cat("miRNA data error:", e$message, "\n\n"))

# 3. Protein Data
tryCatch({
    load(file.path(data_dir, "protein_data.RData"))
    cat("=== PROTEIN DATA ===\n")
    cat("Sample names (first 5):", paste(head(colnames(protein_matrix), 5), collapse=", "), "\n")
    cat("Patient IDs from substr(1,12):", paste(head(unique(substr(colnames(protein_matrix), 1, 12)), 5), collapse=", "), "\n")
    cat("Total samples:", ncol(protein_matrix), "\n")
    cat("Unique patients:", length(unique(substr(colnames(protein_matrix), 1, 12))), "\n\n")
}, error = function(e) cat("Protein data error:", e$message, "\n\n"))

# 4. Subtype Data
tryCatch({
    load(file.path(data_dir, "subtype_info.RData"))
    cat("=== SUBTYPE DATA ===\n")
    cat("Patient IDs (first 5):", paste(head(subtype_info$patient_id, 5), collapse=", "), "\n")
    cat("Total patients:", nrow(subtype_info), "\n")
    cat("Unique subtypes:", paste(unique(subtype_info$subtype), collapse=", "), "\n\n")
}, error = function(e) cat("Subtype data error:", e$message, "\n\n"))

# 5. Patient ID overlap analysis
tryCatch({
    load(file.path(data_dir, "mrna_data.RData"))
    load(file.path(data_dir, "mirna_data.RData"))
    load(file.path(data_dir, "protein_data.RData"))
    load(file.path(data_dir, "subtype_info.RData"))
    
    cat("=== PATIENT OVERLAP ANALYSIS ===\n")
    
    mrna_patients <- unique(substr(colnames(mrna_counts), 1, 12))
    mirna_patients <- unique(substr(colnames(mirna_counts), 1, 12))
    protein_patients <- unique(substr(colnames(protein_matrix), 1, 12))
    subtype_patients <- unique(subtype_info$patient_id)
    
    cat("mRNA patients:", length(mrna_patients), "\n")
    cat("miRNA patients:", length(mirna_patients), "\n")
    cat("protein patients:", length(protein_patients), "\n")
    cat("subtype patients:", length(subtype_patients), "\n\n")
    
    common_all <- Reduce(intersect, list(mrna_patients, mirna_patients, protein_patients))
    common_with_subtype <- intersect(common_all, subtype_patients)
    
    cat("Common across all three data types:", length(common_all), "\n")
    cat("Common with subtype info:", length(common_with_subtype), "\n")
    
    if (length(common_all) > 0) {
        cat("Example common patients:", paste(head(common_all, 3), collapse=", "), "\n")
    } else {
        cat("NO COMMON PATIENTS FOUND!\n")
        
        # Show format mismatches
        cat("\nFormat analysis:\n")
        cat("mRNA TCGA format:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", mrna_patients)), "/", length(mrna_patients), "\n")
        cat("miRNA TCGA format:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", mirna_patients)), "/", length(mirna_patients), "\n")
        cat("protein TCGA format:", sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", protein_patients)), "/", length(protein_patients), "\n")
        
        # Show first few from each
        cat("\nFirst 5 from each:\n")
        cat("mRNA:", paste(head(mrna_patients, 5), collapse=", "), "\n")
        cat("miRNA:", paste(head(mirna_patients, 5), collapse=", "), "\n")
        cat("protein:", paste(head(protein_patients, 5), collapse=", "), "\n")
        cat("subtype:", paste(head(subtype_patients, 5), collapse=", "), "\n")
    }
    
}, error = function(e) cat("Overlap analysis error:", e$message, "\n"))