# Multi-omics Data Preprocessing Script
# Standardizes and prepares multi-omics data for DIABLO analysis

suppressMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(edgeR)
  library(mixOmics)
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

# Function to match samples across data types
match_samples <- function(data_list, subtype_info) {
  cat("Matching samples across data types...\n")
  
  # Get common samples
  sample_names <- colnames(data_list[[1]])
  for (i in 2:length(data_list)) {
    if (!is.null(data_list[[i]])) {
      sample_names <- intersect(sample_names, colnames(data_list[[i]]))
    }
  }
  
  # Filter all data types for common samples
  matched_data <- list()
  for (i in 1:length(data_list)) {
    if (!is.null(data_list[[i]])) {
      matched_data[[i]] <- data_list[[i]][, sample_names]
    } else {
      matched_data[[i]] <- NULL
    }
  }
  names(matched_data) <- names(data_list)
  
  # Filter subtype info for common samples
  matched_subtypes <- subtype_info[subtype_info$sample_id %in% sample_names, ]
  
  cat("Matched", length(sample_names), "samples across data types\n")
  return(list(data = matched_data, subtypes = matched_subtypes))
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

# Process test data
cat("Processing test data...\n")
test_processed <- preprocess_multiomics_data(test_data, "test")

# Save processed data
save(train_processed, file = file.path(multiomics_dir, "train_processed.RData"))
save(test_processed, file = file.path(multiomics_dir, "test_processed.RData"))

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