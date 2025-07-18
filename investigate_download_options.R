#!/usr/bin/env Rscript

# Investigation script to explore download options and metadata mapping for miRNA data
# This will examine different approaches to get TCGA-formatted sample IDs

suppressMessages({
  library(TCGAbiolinks)
  library(dplyr)
})

# Create logs directory
dir.create("logs", showWarnings = FALSE)

# Investigation log file
log_file <- "logs/mirna_download_investigation.txt"
sink(log_file)

cat("=== miRNA Download Options Investigation ===\n")
cat("Timestamp:", Sys.time(), "\n\n")

# =============================================================================
# STEP 1: Explore different miRNA download options
# =============================================================================

cat("1. Exploring different miRNA download approaches...\n")

# Option A: Current approach - individual files approach
cat("\n--- Option A: Current individual files approach ---\n")
tryCatch({
  query_mirna_individual <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification"
  )
  
  cat("Individual files approach - files found:", length(query_mirna_individual$results[[1]]$id), "\n")
  
  # Check if cases are available
  if ("cases" %in% names(query_mirna_individual$results[[1]])) {
    cat("Cases metadata available: YES\n")
    
    # Extract first few cases to examine structure
    cases_data <- query_mirna_individual$results[[1]]$cases
    if (length(cases_data) > 0) {
      cat("First case structure:\n")
      str(cases_data[[1]][[1]])
    }
  } else {
    cat("Cases metadata available: NO\n")
  }
  
}, error = function(e) {
  cat("Error in individual files approach:", e$message, "\n")
})

# Option B: Using GDCprepare with proper metadata extraction
cat("\n--- Option B: GDCprepare approach ---\n")
tryCatch({
  query_mirna_prepare <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling"  # Try specific workflow
  )
  
  cat("GDCprepare approach - files found:", length(query_mirna_prepare$results[[1]]$id), "\n")
  
  # Check metadata structure
  files_info <- query_mirna_prepare$results[[1]]
  available_cols <- names(files_info)
  cat("Available metadata columns:", paste(available_cols, collapse = ", "), "\n")
  
  # Look for specific columns that might contain TCGA IDs
  tcga_columns <- c("submitter_id", "sample_id", "case_id", "barcode", "aliquot_id")
  found_tcga_cols <- tcga_columns[tcga_columns %in% available_cols]
  cat("TCGA-relevant columns found:", paste(found_tcga_cols, collapse = ", "), "\n")
  
}, error = function(e) {
  cat("Error in GDCprepare approach:", e$message, "\n")
})

# =============================================================================
# STEP 2: Examine manifest-based approach
# =============================================================================

cat("\n2. Exploring manifest-based approach...\n")

tryCatch({
  # Create a manifest query
  query_manifest <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification"
  )
  
  # Get the manifest
  manifest <- GDCquery_Maf(query_manifest, save = FALSE)
  
  if (!is.null(manifest) && nrow(manifest) > 0) {
    cat("Manifest approach - entries found:", nrow(manifest), "\n")
    cat("Manifest columns:", paste(colnames(manifest), collapse = ", "), "\n")
    
    # Look for TCGA sample IDs in manifest
    tcga_cols_manifest <- c("submitter_id", "sample_id", "case_id", "barcode")
    found_in_manifest <- tcga_cols_manifest[tcga_cols_manifest %in% colnames(manifest)]
    
    if (length(found_in_manifest) > 0) {
      cat("TCGA columns in manifest:", paste(found_in_manifest, collapse = ", "), "\n")
      
      # Show examples
      cat("\nFirst 3 manifest entries:\n")
      print(head(manifest[, c("id", found_in_manifest)], 3))
    }
  }
  
}, error = function(e) {
  cat("Error in manifest approach:", e$message, "\n")
})

# =============================================================================
# STEP 3: Alternative query approaches
# =============================================================================

cat("\n3. Testing alternative query parameters...\n")

# Approach 1: Different workflow type
tryCatch({
  query_workflow <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling"
  )
  
  if (length(query_workflow$results[[1]]$id) > 0) {
    cat("BCGSC workflow - files:", length(query_workflow$results[[1]]$id), "\n")
    
    # Check metadata structure
    files_info <- query_workflow$results[[1]]
    if ("cases" %in% names(files_info)) {
      case_example <- files_info$cases[[1]][[1]]
      cat("Sample ID from case:")
      if ("samples" %in% names(case_example) && length(case_example$samples) > 0) {
        cat(case_example$samples[[1]]$submitter_id, "\n")
      }
    }
  }
  
}, error = function(e) {
  cat("BCGSC workflow error:", e$message, "\n")
})

# Approach 2: Using legacy data
tryCatch({
  query_legacy <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    legacy = TRUE
  )
  
  cat("Legacy approach - files:", length(query_legacy$results[[1]]$id), "\n")
  
}, error = function(e) {
  cat("Legacy approach error:", e$message, "\n")
})

# =============================================================================
# STEP 4: Direct API approach to get mapping
# =============================================================================

cat("\n4. Testing direct API approach for mapping...\n")

tryCatch({
  # Get basic query
  query_basic <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification"
  )
  
  # Extract file IDs and attempt to get metadata
  file_ids <- query_basic$results[[1]]$id[1:min(5, length(query_basic$results[[1]]$id))]
  
  cat("Examining first", length(file_ids), "file(s) for metadata:\n")
  
  for (file_id in file_ids) {
    cat("\nFile ID:", file_id, "\n")
    
    # Try to get associated case/sample information
    tryCatch({
      # Use GDC API to get metadata
      metadata_url <- paste0("https://api.gdc.cancer.gov/files/", file_id)
      
      # This would require httr package, so we'll simulate the approach
      # For now, examine what's available in the query results
      
      files_info <- query_basic$results[[1]]
      row_idx <- which(files_info$id == file_id)
      
      if (length(row_idx) > 0) {
        # Check all available columns for this file
        cat("  Available metadata:\n")
        for (col_name in names(files_info)) {
          if (!is.null(files_info[[col_name]][row_idx])) {
            cat("    ", col_name, ":", 
                if(is.character(files_info[[col_name]][row_idx])) 
                  substr(as.character(files_info[[col_name]][row_idx]), 1, 50)
                else 
                  "[non-character]", "\n")
          }
        }
      }
      
    }, error = function(e) {
      cat("  Error getting metadata:", e$message, "\n")
    })
  }
  
}, error = function(e) {
  cat("API approach error:", e$message, "\n")
})

# =============================================================================
# STEP 5: Summary and recommendations
# =============================================================================

cat("\n5. Summary of investigation:\n")
cat("\nCurrent situation:")
cat("- miRNA data has UUID format sample IDs (1207 samples)")
cat("- mRNA data has TCGA format sample IDs (1231 samples)")
cat("- No overlap between datasets due to format mismatch")

cat("\nPotential solutions:")
cat("1. Use GDCprepare() with proper metadata extraction")
cat("2. Create UUID-to-TCGA mapping table")
cat("3. Use manifest-based download with mapping")
cat("4. Post-process UUIDs to extract TCGA barcodes")

cat("\nNext steps needed:")
cat("- Implement proper metadata extraction in miRNA download")
cat("- Create UUID-to-TCGA mapping function")
cat("- Test mapping with sample data")

sink()

cat("Investigation complete! Results saved to:", log_file, "\n")