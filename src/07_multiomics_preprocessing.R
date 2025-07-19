# Multi-omics Data Preprocessing Script
# Standardizes and prepares multi-omics data for DIABLO analysis

suppressMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(edgeR)
  library(mixOmics)
})

# Ensure logs directory exists and setup safe output capture
if (!dir.exists("logs")) {
  dir.create("logs", showWarnings = FALSE)
}

tryCatch({
  # Check if sink is already active
  if (sink.number() == 0) {
    sink("logs/tcga_pipeline_output.log", append = TRUE)
    sink("logs/tcga_pipeline_output.log", type = "message", append = TRUE)
  }
}, error = function(e) {
  warning("Could not set up output capture: ", e$message)
})

cat("Loading multi-omics preprocessing functions...\n")

# Load the downloaded multi-omics data
multiomics_dir <- file.path(data_dir, "multiomics")
load(file.path(multiomics_dir, "train_data.RData"))
load(file.path(multiomics_dir, "test_data.RData"))

# Function to convert ENSEMBL IDs to gene symbols for mRNA data
convert_ensembl_to_symbols <- function(count_matrix) {
  cat("Converting ENSEMBL IDs to gene symbols...\n")
  
  # Remove version numbers from ENSEMBL IDs
  ensembl_ids <- gsub("\\..*", "", rownames(count_matrix))
  
  # Map ENSEMBL IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db, 
                        keys = ensembl_ids,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  
  # Keep only genes with valid symbols
  valid_genes <- !is.na(gene_symbols)
  gene_symbols <- gene_symbols[valid_genes]
  count_matrix <- count_matrix[valid_genes, ]
  
  # Handle duplicate gene symbols
  duplicated_symbols <- duplicated(gene_symbols)
  gene_symbols <- gene_symbols[!duplicated_symbols]
  count_matrix <- count_matrix[!duplicated_symbols, ]
  
  # Set gene symbols as row names
  rownames(count_matrix) <- gene_symbols
  
  cat("Converted", nrow(count_matrix), "genes with valid symbols\n")
  return(count_matrix)
}

# Function to preprocess mRNA data
preprocess_mrna <- function(mrna_counts) {
  cat("Preprocessing mRNA data...\n")
  
  # Convert ENSEMBL IDs to gene symbols
  mrna_symbols <- convert_ensembl_to_symbols(mrna_counts)
  
  # Filter low-expressed genes
  keep_genes <- rowSums(mrna_symbols > 1) >= 10
  mrna_filtered <- mrna_symbols[keep_genes, ]
  
  # CPM normalization and log transformation
  mrna_cpm <- cpm(mrna_filtered, log = TRUE, prior.count = 1)
  
  # Select most variable genes
  gene_vars <- apply(mrna_cpm, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
  mrna_final <- mrna_cpm[top_genes, ]
  
  cat("mRNA data: ", nrow(mrna_final), "genes x", ncol(mrna_final), "samples\n")
  return(mrna_final)
}

# Function to preprocess miRNA data
preprocess_mirna <- function(mirna_counts) {
  cat("Preprocessing miRNA data...\n")
  
  # Filter low-expressed miRNAs
  keep_mirnas <- rowSums(mirna_counts > 1) >= 10
  mirna_filtered <- mirna_counts[keep_mirnas, ]
  
  # CPM normalization and log transformation
  mirna_cpm <- cpm(mirna_filtered, log = TRUE, prior.count = 1)
  
  # Select most variable miRNAs
  mirna_vars <- apply(mirna_cpm, 1, var)
  top_mirnas <- names(sort(mirna_vars, decreasing = TRUE))[1:200]
  mirna_final <- mirna_cpm[top_mirnas, ]
  
  cat("miRNA data: ", nrow(mirna_final), "miRNAs x", ncol(mirna_final), "samples\n")
  return(mirna_final)
}

# Function to preprocess protein data
preprocess_protein <- function(protein_matrix) {
  cat("Preprocessing protein data...\n")
  
  # Remove proteins with too many missing values
  missing_prop <- apply(protein_matrix, 1, function(x) sum(is.na(x)) / length(x))
  keep_proteins <- missing_prop < 0.3
  protein_filtered <- protein_matrix[keep_proteins, ]
  
  # Impute missing values with median
  for (i in 1:nrow(protein_filtered)) {
    missing_idx <- is.na(protein_filtered[i, ])
    if (sum(missing_idx) > 0) {
      protein_filtered[i, missing_idx] <- median(protein_filtered[i, ], na.rm = TRUE)
    }
  }
  
  # Scale and center the data
  protein_scaled <- t(scale(t(protein_filtered)))
  
  cat("Protein data: ", nrow(protein_scaled), "proteins x", ncol(protein_scaled), "samples\n")
  return(protein_scaled)
}

# Function to harmonize sample names across omics types
harmonize_sample_names <- function(data_list, subtype_info) {
  log_info("Harmonizing sample names across omics types...")
  
  # Extract patient IDs from sample barcodes
  extract_patient_id <- function(sample_barcode) {
    # Extract patient ID (TCGA-XX-XXXX format)
    if (grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_barcode)) {
      return(substr(sample_barcode, 1, 12))
    } else {
      return(sample_barcode)  # fallback
    }
  }
  
  # Create mapping from original sample names to patient IDs
  sample_mapping <- list()
  
  for (block_name in names(data_list)) {
    if (!is.null(data_list[[block_name]])) {
      original_names <- colnames(data_list[[block_name]])
      patient_ids <- sapply(original_names, extract_patient_id)
      sample_mapping[[block_name]] <- data.frame(
        original_name = original_names,
        patient_id = patient_ids,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Find common patient IDs across all omics types
  all_patient_ids <- lapply(sample_mapping, function(x) unique(x$patient_id))
  common_patients <- Reduce(intersect, all_patient_ids)
  
  log_info("Found %d common patients across all omics types", length(common_patients))
  
  # Harmonize data matrices to use patient-level sample names
  harmonized_data <- list()
  for (block_name in names(data_list)) {
    if (!is.null(data_list[[block_name]])) {
      # Map original sample names to patient IDs
      mapping <- sample_mapping[[block_name]]
      
      # Keep only samples from common patients
      keep_samples <- mapping$patient_id %in% common_patients
      filtered_data <- data_list[[block_name]][, keep_samples]
      
      # Rename columns to use patient IDs
      new_names <- mapping$patient_id[keep_samples]
      colnames(filtered_data) <- new_names
      
      harmonized_data[[block_name]] <- filtered_data
      
      log_info("Harmonized %s: %d samples (%d unique patients)", 
               block_name, ncol(filtered_data), length(unique(new_names)))
    }
  }
  
  # Update subtype info to use patient IDs
  harmonized_subtypes <- subtype_info
  harmonized_subtypes$sample_id <- sapply(harmonized_subtypes$sample_id, extract_patient_id)
  
  # Filter to only include common patients
  harmonized_subtypes <- harmonized_subtypes[harmonized_subtypes$sample_id %in% common_patients, ]
  
  # Remove duplicates (keep first occurrence per patient)
  harmonized_subtypes <- harmonized_subtypes[!duplicated(harmonized_subtypes$sample_id), ]
  
  log_info("Harmonized subtypes: %d patients with subtype info", nrow(harmonized_subtypes))
  
  return(list(
    data = harmonized_data,
    subtypes = harmonized_subtypes,
    sample_mapping = sample_mapping
  ))
}

# Function to create design matrix for DIABLO
create_design_matrix <- function(n_blocks) {
  cat("Creating design matrix for DIABLO...\n")
  
  # Create design matrix with 0.1 correlation between all blocks
  design <- matrix(0.1, nrow = n_blocks, ncol = n_blocks)
  diag(design) <- 0
  
  # Set block names
  if (n_blocks == 3) {
    rownames(design) <- colnames(design) <- c("mRNA", "miRNA", "protein")
  } else if (n_blocks == 2) {
    rownames(design) <- colnames(design) <- c("mRNA", "miRNA")
  }
  
  cat("Design matrix created with", n_blocks, "blocks\n")
  return(design)
}


match_samples <- function(data_list, subtype_info) {
  log_info("Matching samples across data types...")
  
  # First harmonize sample names to use patient-level identifiers
  harmonized_result <- harmonize_sample_names(data_list, subtype_info)
  
  # Verify that all data types have the same sample names
  sample_names <- NULL
  for (block_name in names(harmonized_result$data)) {
    if (!is.null(harmonized_result$data[[block_name]])) {
      current_names <- colnames(harmonized_result$data[[block_name]])
      if (is.null(sample_names)) {
        sample_names <- current_names
      } else if (!identical(sample_names, current_names)) {
        log_error("Sample name mismatch detected after harmonization: %s vs %s", 
                 head(sample_names, 3), head(current_names, 3))
        stop("Sample name harmonization failed")
      }
    }
  }
  
  # Verify subtype alignment
  subtype_samples <- harmonized_result$subtypes$sample_id
  if (!identical(sort(sample_names), sort(subtype_samples))) {
    log_error("Subtype samples do not match harmonized data samples")
    log_error("Data samples: %s", paste(head(sample_names, 5), collapse=", "))
    log_error("Subtype samples: %s", paste(head(subtype_samples, 5), collapse=", "))
    stop("Sample alignment failed")
  }
  
  log_info("Successfully matched %d samples across all data types", length(sample_names))
  return(list(data = harmonized_result$data, subtypes = harmonized_result$subtypes))
}

# Main preprocessing function
preprocess_multiomics_data <- function(raw_data, dataset_name) {
  cat("=== Preprocessing", dataset_name, "data ===\n")
  
  # Preprocess each data type
  processed_data <- list()
  
  # mRNA preprocessing
  processed_data$mRNA <- preprocess_mrna(raw_data$mrna)
  
  # miRNA preprocessing
  processed_data$miRNA <- preprocess_mirna(raw_data$mirna)
  
  # Protein preprocessing (if available)
  if (!is.null(raw_data$protein)) {
    processed_data$protein <- preprocess_protein(raw_data$protein)
  } else {
    processed_data$protein <- NULL
  }
  
  # Match samples across data types
  matched_result <- match_samples(processed_data, raw_data$subtypes)
  
  # Create design matrix
  n_blocks <- sum(!sapply(matched_result$data, is.null))
  design_matrix <- create_design_matrix(n_blocks)
  
  # Prepare final data for mixOmics
  final_data <- list()
  for (block_name in names(matched_result$data)) {
    if (!is.null(matched_result$data[[block_name]])) {
      # Transpose for mixOmics (samples as rows, features as columns)
      final_data[[block_name]] <- t(matched_result$data[[block_name]])
    }
  }
  
  # Create outcome vector
  Y <- as.factor(matched_result$subtypes$subtype)
  names(Y) <- matched_result$subtypes$sample_id
  
  # DEBUG: Log data before transposing (with harmonized names)
  log_info("=== DEBUG: Before Transposing (Harmonized) ===")
  for (block_name in names(matched_result$data)) {
    if (!is.null(matched_result$data[[block_name]])) {
      log_info("%s data (features x samples): dimensions %dx%d", 
               block_name, 
               nrow(matched_result$data[[block_name]]), 
               ncol(matched_result$data[[block_name]]))
      log_info("First 5 rows x 5 cols:\n%s", 
               paste(capture.output(head(matched_result$data[[block_name]][1:5, 1:5])), collapse="\n"))
      log_info("Column names (samples, first 10): %s", 
               paste(head(colnames(matched_result$data[[block_name]]), 10), collapse=", "))
      log_info("Row names (features, first 10): %s", 
               paste(head(rownames(matched_result$data[[block_name]]), 10), collapse=", "))
    }
  }

  # Prepare final data for mixOmics
  final_data <- list()
  for (block_name in names(matched_result$data)) {
    if (!is.null(matched_result$data[[block_name]])) {
      # Transpose for mixOmics (samples as rows, features as columns)
      final_data[[block_name]] <- t(matched_result$data[[block_name]])
      
      # DEBUG: Log data after transposing (with harmonized names)
      log_info("=== DEBUG: After Transposing %s (Harmonized) ===", block_name)
      log_info("Dimensions: %dx%d", nrow(final_data[[block_name]]), ncol(final_data[[block_name]]))
      log_info("First 5 rows x 5 cols:\n%s", 
               paste(capture.output(head(final_data[[block_name]][1:5, 1:5])), collapse="\n"))
      log_info("Row names (samples, first 10): %s", 
               paste(head(rownames(final_data[[block_name]]), 10), collapse=", "))
      log_info("Column names (features, first 10): %s", 
               paste(head(colnames(final_data[[block_name]]), 10), collapse=", "))
    }
  }

  # Create outcome vector
  Y <- as.factor(matched_result$subtypes$subtype)
  names(Y) <- matched_result$subtypes$sample_id
  
  # DEBUG: Log outcome vector information
  log_info("=== DEBUG: Outcome Vector (Harmonized) ===")
  log_info("Y vector length: %d", length(Y))
  log_info("Y names (samples, first 10): %s", 
           paste(head(names(Y), 10), collapse=", "))
  log_info("Y levels: %s", paste(levels(Y), collapse=", "))
  log_info("Y distribution:\n%s", paste(capture.output(table(Y)), collapse="\n"))
  
  # DEBUG: Check row name alignment across all blocks
  log_info("=== DEBUG: Sample Alignment Check (Harmonized) ===")
  for (block_name in names(final_data)) {
    if (!is.null(final_data[[block_name]])) {
      match_status <- all(rownames(final_data[[block_name]]) == names(Y))
      log_info("%s row names match Y names: %s", block_name, match_status)
      if (!match_status) {
        log_error("Mismatch detected! Sample names do not align")
        log_error("Data samples: %s", paste(head(rownames(final_data[[block_name]]), 5), collapse=", "))
        log_error("Y samples: %s", paste(head(names(Y), 5), collapse=", "))
        stop("Sample alignment verification failed")
      }
    }
  }

  return(list(
    X = final_data,
    Y = Y,
    design = design_matrix,
    subtypes = matched_result$subtypes
  ))
}

# Process training data
cat("Processing training data...\n")
train_processed <- preprocess_multiomics_data(train_data, "training")

# DEBUG: Print training processed data structure
save(train_processed, file = file.path(multiomics_dir, "train_processed.RData"))
test_processed <- preprocess_multiomics_data(test_data, "test")

# DEBUG: Log processed data structures
log_info("=== DEBUG: Training Processed Data ===")
for (block_name in names(train_processed$X)) {
  if (!is.null(train_processed$X[[block_name]])) {
    log_data_summary(train_processed$X[[block_name]], paste(block_name, "training data"), "matrix")
    log_info("First 5 rows x 5 cols:\n%s", 
             paste(capture.output(head(train_processed$X[[block_name]][1:5, 1:5])), collapse="\n"))
    log_info("Row names (samples, first 10): %s", 
             paste(head(rownames(train_processed$X[[block_name]]), 10), collapse=", "))
  }
}

log_info("=== DEBUG: Test Processed Data ===")
for (block_name in names(test_processed$X)) {
  if (!is.null(test_processed$X[[block_name]])) {
    log_data_summary(test_processed$X[[block_name]], paste(block_name, "test data"), "matrix")
    log_info("First 5 rows x 5 cols:\n%s", 
             paste(capture.output(head(test_processed$X[[block_name]][1:5, 1:5])), collapse="\n"))
    log_info("Row names (samples, first 10): %s", 
             paste(head(rownames(test_processed$X[[block_name]]), 10), collapse=", "))
  }
}

# DEBUG: Cross-validation of sample alignment
log_info("=== DEBUG: Final Sample Alignment (Harmonized) ===")
log_info("Training Y names vs X row names:")
for (block_name in names(train_processed$X)) {
  if (!is.null(train_processed$X[[block_name]])) {
    match_status <- all(rownames(train_processed$X[[block_name]]) == names(train_processed$Y))
    log_info("  %s: %s", block_name, match_status)
    if (!match_status) {
      log_error("Training set alignment issue detected in %s", block_name)
      log_error("  Data samples: %s", paste(head(rownames(train_processed$X[[block_name]]), 3), collapse=", "))
      log_error("  Y samples: %s", paste(head(names(train_processed$Y), 3), collapse=", "))
      stop("Training set alignment verification failed")
    }
  }
}

log_info("Test Y names vs X row names:")
for (block_name in names(test_processed$X)) {
  if (!is.null(test_processed$X[[block_name]])) {
    match_status <- all(rownames(test_processed$X[[block_name]]) == names(test_processed$Y))
    log_info("  %s: %s", block_name, match_status)
    if (!match_status) {
      log_error("Test set alignment issue detected in %s", block_name)
      log_error("  Data samples: %s", paste(head(rownames(test_processed$X[[block_name]]), 3), collapse=", "))
      log_error("  Y samples: %s", paste(head(names(test_processed$Y), 3), collapse=", "))
      stop("Test set alignment verification failed")
    }
  }
}

log_info("=== Sample Name Harmonization Summary ===")
log_info("All omics blocks now use consistent patient-level sample names")
log_info("Block sPLS-DA requirements for sample name alignment should be satisfied")

# Create summary statistics
create_preprocessing_summary <- function(processed_data, dataset_name) {
  summary_stats <- list()
  
  for (block_name in names(processed_data$X)) {
    if (!is.null(processed_data$X[[block_name]])) {
      summary_stats[[block_name]] <- list(
        n_samples = nrow(processed_data$X[[block_name]]),
        n_features = ncol(processed_data$X[[block_name]]),
        feature_range = range(processed_data$X[[block_name]]),
        mean_expression = mean(processed_data$X[[block_name]]),
        sd_expression = sd(processed_data$X[[block_name]])
      )
    }
  }
  
  summary_stats$outcome <- table(processed_data$Y)
  summary_stats$dataset <- dataset_name
  
  return(summary_stats)
}

# Generate summaries
train_summary <- create_preprocessing_summary(train_processed, "training")
test_summary <- create_preprocessing_summary(test_processed, "test")

# Save summaries
save(train_summary, test_summary, file = file.path(multiomics_dir, "preprocessing_summary.RData"))

# Print summary
cat("\n=== Preprocessing Summary ===\n")
cat("Training data:\n")
for (block in names(train_summary)) {
  if (block != "outcome" && block != "dataset") {
    cat("  ", block, ": ", train_summary[[block]]$n_features, "features x", 
        train_summary[[block]]$n_samples, "samples\n")
  }
}
cat("  Outcome distribution:", paste(names(train_summary$outcome), 
                                   train_summary$outcome, sep = "=", collapse = ", "), "\n")

cat("Test data:\n")
for (block in names(test_summary)) {
  if (block != "outcome" && block != "dataset") {
    cat("  ", block, ": ", test_summary[[block]]$n_features, "features x", 
        test_summary[[block]]$n_samples, "samples\n")
  }
}
cat("  Outcome distribution:", paste(names(test_summary$outcome), 
                                   test_summary$outcome, sep = "=", collapse = ", "), "\n")

cat("Multi-omics preprocessing completed successfully!\n")

# Close output capture - only close if we opened it
tryCatch({
  if (sink.number("output") > 0) {
    sink(type = "output")
  }
  if (sink.number("message") > 0) {
    sink(type = "message")
  }
}, error = function(e) {
  warning("Could not close output capture: ", e$message)
})