# Multi-omics Data Download Script
# Downloads TCGA breast cancer data (mRNA, miRNA, protein) following DIABLO methodology

suppressMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(SummarizedExperiment)
})

# Load logging utilities
source(file.path("src", "utils", "logging.R"))

# Initialize logging if not already done
if (is.null(.log_config)) {
  initialize_logging()
}

log_info("=== Multi-omics Data Download Module ===")
log_info("Loading multi-omics data download functions...")

# Create multi-omics data directory
multiomics_dir <- file.path(data_dir, "multiomics")
dir.create(multiomics_dir, recursive = TRUE, showWarnings = FALSE)

# Download mRNA data for breast cancer
download_mrna_data <- function() {
  operation_id <- log_operation_start("mRNA Data Download", "TCGA-BRCA Gene Expression Quantification")
  
  # Check if processed data already exists
  save_path <- file.path(multiomics_dir, "mrna_data.RData")
  if (file.exists(save_path)) {
    log_info("Found existing mRNA data, loading from file...")
    load(save_path)
    log_data_summary(mrna_counts, "mRNA counts matrix (loaded)")
    log_data_summary(mrna_samples, "mRNA sample metadata (loaded)", "data.frame")
    log_operation_end(operation_id, "Success (loaded existing data)")
    return(list(counts = mrna_counts, samples = mrna_samples))
  }
  
  result <- with_error_handling({
    # Query mRNA data
    log_info("Querying mRNA data from GDC...")
    query_mrna <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    log_info("Found %d mRNA files for download", length(query_mrna$results[[1]]$id))
    
    # Download the data
    log_info("Downloading mRNA data files...")
    GDCdownload(query_mrna, method = "api", files.per.chunk = 10)
    
    # Prepare the data
    log_info("Preparing mRNA data...")
    mrna_data <- GDCprepare(query_mrna)
    
    # Extract counts matrix
    mrna_counts <- assay(mrna_data)
    log_data_summary(mrna_counts, "mRNA counts matrix")
    
    # Extract sample information
    mrna_samples <- as.data.frame(colData(mrna_data))
    log_data_summary(mrna_samples, "mRNA sample metadata", "data.frame")
    
    # Save data
    save(mrna_counts, mrna_samples, file = save_path)
    log_file_operation("saved", save_path, TRUE, "mRNA data")
    
    # Clean up temporary files
    cleanup_temp_files()
    
    list(counts = mrna_counts, samples = mrna_samples)
  }, "mRNA Data Download")
  
  log_operation_end(operation_id, if(!is.null(result)) "Success" else "Failed")
  return(result)
}

# Helper function to process individual miRNA files
process_mirna_files <- function(query_mirna) {
  log_info("Processing individual miRNA quantification files...")
  
  # Get the file information from the query
  files_info <- query_mirna$results[[1]]
  
  # Find all downloaded miRNA files
  mirna_dir <- file.path("GDCdata", "TCGA-BRCA", "Transcriptome_Profiling", "miRNA_Expression_Quantification")
  
  # Initialize data structures
  mirna_counts_list <- list()
  mirna_ids <- NULL
  
  # Process each sample directory
  for (i in 1:nrow(files_info)) {
    file_id <- files_info$id[i]
    sample_dir <- file.path(mirna_dir, file_id)
    
    if (dir.exists(sample_dir)) {
      # Find the quantification file in the directory
      quant_files <- list.files(sample_dir, pattern = "\\.mirbase21\\.mirnas\\.quantification\\.txt$", full.names = TRUE)
      
      if (length(quant_files) > 0) {
        quant_file <- quant_files[1]
        
        # Read the quantification file
        mirna_data <- read.delim(quant_file, header = TRUE, stringsAsFactors = FALSE)
        
        # Extract read counts (second column)
        counts <- mirna_data$read_count
        names(counts) <- mirna_data$miRNA_ID
        
        # Store in list with file UUID as name
        mirna_counts_list[[file_id]] <- counts
        
        # Set miRNA IDs if not already set
        if (is.null(mirna_ids)) {
          mirna_ids <- mirna_data$miRNA_ID
        }
      }
    }
    
    # Log progress periodically
    if (i %% 100 == 0) {
      log_progress(i, nrow(files_info), "Processing miRNA files")
    }
  }
  
  log_info("Successfully processed %d miRNA files", length(mirna_counts_list))
  
  # Create the count matrix
  log_info("Building miRNA count matrix...")
  sample_ids <- names(mirna_counts_list)
  
  # Initialize matrix
  mirna_counts <- matrix(0, nrow = length(mirna_ids), ncol = length(sample_ids))
  rownames(mirna_counts) <- mirna_ids
  colnames(mirna_counts) <- sample_ids
  
  # Fill the matrix
  for (sample_id in sample_ids) {
    sample_counts <- mirna_counts_list[[sample_id]]
    # Match miRNA IDs and fill counts
    matched_ids <- match(mirna_ids, names(sample_counts))
    mirna_counts[, sample_id] <- ifelse(is.na(matched_ids), 0, sample_counts[matched_ids])
  }
  
  # Try to map UUIDs to TCGA sample IDs
  log_info("Mapping file UUIDs to TCGA sample IDs...")
  
  # Use the file information to map UUIDs to sample IDs
  if ("cases" %in% names(files_info)) {
    # Extract sample IDs from the cases information
    sample_mapping <- data.frame(
      file_id = files_info$id,
      sample_id = sapply(files_info$cases, function(x) {
        if (length(x) > 0 && "samples" %in% names(x[[1]])) {
          x[[1]]$samples[[1]]$submitter_id
        } else {
          NA
        }
      }),
      stringsAsFactors = FALSE
    )
    
    # Update column names where mapping is available
    for (i in 1:nrow(sample_mapping)) {
      if (!is.na(sample_mapping$sample_id[i])) {
        col_idx <- which(colnames(mirna_counts) == sample_mapping$file_id[i])
        if (length(col_idx) > 0) {
          colnames(mirna_counts)[col_idx] <- sample_mapping$sample_id[i]
        }
      }
    }
    
    # Create sample metadata
    mirna_samples <- data.frame(
      submitter_id = colnames(mirna_counts),
      sample_type = "Primary Tumor",
      stringsAsFactors = FALSE
    )
    rownames(mirna_samples) <- colnames(mirna_counts)
    
  } else {
    # Fallback: use file IDs as sample IDs
    log_warn("Could not map file UUIDs to sample IDs, using file IDs")
    mirna_samples <- data.frame(
      submitter_id = colnames(mirna_counts),
      sample_type = "Primary Tumor",
      stringsAsFactors = FALSE
    )
    rownames(mirna_samples) <- colnames(mirna_counts)
  }
  
  return(list(counts = mirna_counts, samples = mirna_samples))
}

# Download miRNA data for breast cancer
download_mirna_data <- function() {
  operation_id <- log_operation_start("miRNA Data Download", "TCGA-BRCA miRNA Expression Quantification")
  
  # Check if processed data already exists
  save_path <- file.path(multiomics_dir, "mirna_data.RData")
  if (file.exists(save_path)) {
    log_info("Found existing miRNA data, loading from file...")
    load(save_path)
    log_data_summary(mirna_counts, "miRNA counts matrix (loaded)")
    log_data_summary(mirna_samples, "miRNA sample metadata (loaded)", "data.frame")
    log_operation_end(operation_id, "Success (loaded existing data)")
    return(list(counts = mirna_counts, samples = mirna_samples))
  }
  
  result <- with_error_handling({
    # Query miRNA data
    log_info("Querying miRNA data from GDC...")
    query_mirna <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "miRNA Expression Quantification"
    )
    
    log_info("Found %d miRNA files for download", length(query_mirna$results[[1]]$id))
    
    # Download the data
    log_info("Downloading miRNA data files...")
    GDCdownload(query_mirna, method = "api", files.per.chunk = 10)
    
    # Process the individual miRNA files
    log_info("Processing miRNA files...")
    mirna_result <- process_mirna_files(query_mirna)
    mirna_counts <- mirna_result$counts
    mirna_samples <- mirna_result$samples
    
    log_data_summary(mirna_counts, "miRNA counts matrix")
    log_data_summary(mirna_samples, "miRNA sample metadata", "data.frame")
    
    # Save data
    save(mirna_counts, mirna_samples, file = save_path)
    log_file_operation("saved", save_path, TRUE, "miRNA data")
    
    # Clean up temporary files
    cleanup_temp_files()
    
    list(counts = mirna_counts, samples = mirna_samples)
  }, "miRNA Data Download")
  
  log_operation_end(operation_id, if(!is.null(result)) "Success" else "Failed")
  return(result)
}

# Helper function to clean up temporary download files
cleanup_temp_files <- function() {
  log_info("Cleaning up temporary download files...")
  
  # Pattern for temporary download files (UUID format)
  temp_pattern <- "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
  
  # Find all files matching the pattern
  all_files <- list.files(".", full.names = FALSE)
  temp_files <- all_files[grepl(temp_pattern, all_files)]
  
  if (length(temp_files) > 0) {
    log_info("Found %d temporary files to clean up", length(temp_files))
    
    for (temp_file in temp_files) {
      if (file.exists(temp_file)) {
        # Check if it's a directory or file
        if (file.info(temp_file)$isdir) {
          unlink(temp_file, recursive = TRUE)
        } else {
          file.remove(temp_file)
        }
      }
    }
    
    log_info("Cleaned up %d temporary files", length(temp_files))
  } else {
    log_info("No temporary files found to clean up")
  }
}

# Download protein data for breast cancer
download_protein_data <- function() {
  operation_id <- log_operation_start("Protein Data Download", "TCGA-BRCA Protein Expression Quantification")
  
  # Check if processed data already exists
  save_path <- file.path(multiomics_dir, "protein_data.RData")
  if (file.exists(save_path)) {
    log_info("Found existing protein data, loading from file...")
    load(save_path)
    log_data_summary(protein_matrix, "protein matrix (loaded)")
    log_data_summary(protein_samples, "protein sample metadata (loaded)", "data.frame")
    log_operation_end(operation_id, "Success (loaded existing data)")
    return(list(matrix = protein_matrix, samples = protein_samples))
  }
  
  result <- with_error_handling({
    # Query protein data
    log_info("Querying protein data from GDC...")
    query_protein <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Proteome Profiling",
      data.type = "Protein Expression Quantification"
    )
    
    log_info("Found %d protein files for download", length(query_protein$results[[1]]$id))
    
    # Download the data
    log_info("Downloading protein data files...")
    GDCdownload(query_protein, method = "api", files.per.chunk = 10)
    
    # Prepare the data
    log_info("Preparing protein data...")
    protein_data <- GDCprepare(query_protein)
    
    # Check if protein_data is a SummarizedExperiment or data.frame
    if (is(protein_data, "SummarizedExperiment")) {
      log_info("Protein data returned as SummarizedExperiment, using assay() method")
      # Extract protein matrix using assay() for SummarizedExperiment
      protein_matrix <- assay(protein_data)
      # Extract sample information
      protein_samples <- as.data.frame(colData(protein_data))
    } else if (is(protein_data, "data.frame")) {
      # Handle data.frame format - protein data often comes as data.frame
      log_info("Protein data returned as data.frame, processing accordingly...")
      
      # For protein data.frame, samples are usually columns and proteins are rows
      # Find sample columns (typically start with TCGA- or are numeric/sample identifiers)
      
      # Look for column patterns that might be sample IDs
      sample_cols <- grep("^TCGA-", colnames(protein_data), value = TRUE)
      
      if (length(sample_cols) == 0) {
        # Try alternative patterns - numeric columns or columns with sample-like names
        numeric_cols <- sapply(protein_data, is.numeric)
        sample_cols <- names(numeric_cols)[numeric_cols]
        
        # Remove obvious metadata columns
        metadata_patterns <- c("gene", "protein", "symbol", "id", "name", "description")
        sample_cols <- sample_cols[!grepl(paste(metadata_patterns, collapse = "|"), 
                                         sample_cols, ignore.case = TRUE)]
      }
      
      log_info("Found %d potential sample columns in protein data", length(sample_cols))
      
      if (length(sample_cols) > 0) {
        # Extract protein matrix
        protein_matrix <- as.matrix(protein_data[, sample_cols])
        
        # Set row names - look for gene/protein identifier columns
        if ("Gene" %in% colnames(protein_data)) {
          rownames(protein_matrix) <- protein_data$Gene
        } else if ("gene_name" %in% colnames(protein_data)) {
          rownames(protein_matrix) <- protein_data$gene_name
        } else if ("protein_id" %in% colnames(protein_data)) {
          rownames(protein_matrix) <- protein_data$protein_id
        } else {
          # Use row numbers if no identifier column found
          rownames(protein_matrix) <- paste0("protein_", 1:nrow(protein_matrix))
        }
        
        # Create sample information data.frame
        protein_samples <- data.frame(
          submitter_id = colnames(protein_matrix),
          sample_type = "Primary Tumor",  # Default for TCGA-BRCA
          stringsAsFactors = FALSE
        )
        rownames(protein_samples) <- colnames(protein_matrix)
        
      } else {
        stop("Could not identify sample columns in protein data")
      }
      
    } else {
      stop("Unexpected data format for protein data: ", class(protein_data))
    }
    
    log_data_summary(protein_matrix, "protein matrix")
    log_data_summary(protein_samples, "protein sample metadata", "data.frame")
    
    # Save data
    save(protein_matrix, protein_samples, file = save_path)
    log_file_operation("saved", save_path, TRUE, "protein data")
    
    # Clean up temporary files
    cleanup_temp_files()
    
    list(matrix = protein_matrix, samples = protein_samples)
  }, "Protein Data Download")
  
  log_operation_end(operation_id, if(!is.null(result)) "Success" else "Failed")
  return(result)
}

# Get breast cancer subtype information
get_subtype_info <- function() {
  cat("Downloading breast cancer subtype information...\n")
  
  # Query clinical data
  query_clinical <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR XML"
  )
  
  # Download clinical data
  GDCdownload(query_clinical)
  
  # Prepare clinical data
  clinical_data <- GDCprepare(query_clinical)
  
  # Extract subtype information
  subtype_info <- clinical_data %>%
    select(submitter_id, paper_BRCA_Subtype_PAM50) %>%
    filter(!is.na(paper_BRCA_Subtype_PAM50)) %>%
    rename(
      sample_id = submitter_id,
      subtype = paper_BRCA_Subtype_PAM50
    )
  
  # Save subtype information
  save(subtype_info, file = file.path(multiomics_dir, "subtype_info.RData"))
  cat("Subtype information downloaded and saved.\n")
  
  return(subtype_info)
}

# Create training and test sets following DIABLO methodology
create_train_test_sets <- function(mrna_data, mirna_data, protein_data, subtype_info) {
  cat("Creating training and test sets...\n")
  
  # Get common samples across all data types
  common_samples <- intersect(
    intersect(colnames(mrna_data$counts), colnames(mirna_data$counts)),
    colnames(protein_data$matrix)
  )
  
  # Filter subtype info for common samples
  available_subtypes <- subtype_info[subtype_info$sample_id %in% common_samples, ]
  
  # Focus on main subtypes: Basal, Her2, LumA
  main_subtypes <- available_subtypes[available_subtypes$subtype %in% c("Basal", "Her2", "LumA"), ]
  
  # Create training set (150 samples with all three data types)
  set.seed(42)  # For reproducibility
  train_samples <- main_subtypes %>%
    group_by(subtype) %>%
    sample_n(50) %>%  # 50 samples per subtype
    ungroup() %>%
    pull(sample_id)
  
  # Create test set (70 samples with only mRNA and miRNA)
  remaining_samples <- main_subtypes[!main_subtypes$sample_id %in% train_samples, ]
  test_samples <- remaining_samples %>%
    group_by(subtype) %>%
    sample_n(min(25, n())) %>%  # Up to 25 samples per subtype
    ungroup() %>%
    pull(sample_id)
  
  # Create training data
  train_data <- list(
    mrna = mrna_data$counts[, train_samples],
    mirna = mirna_data$counts[, train_samples],
    protein = protein_data$matrix[, train_samples],
    subtypes = main_subtypes[main_subtypes$sample_id %in% train_samples, ]
  )
  
  # Create test data (missing protein data)
  test_data <- list(
    mrna = mrna_data$counts[, test_samples],
    mirna = mirna_data$counts[, test_samples],
    protein = NULL,  # Missing protein data as in the blog post
    subtypes = main_subtypes[main_subtypes$sample_id %in% test_samples, ]
  )
  
  # Save training and test sets
  save(train_data, file = file.path(multiomics_dir, "train_data.RData"))
  save(test_data, file = file.path(multiomics_dir, "test_data.RData"))
  
  cat("Training set created with", length(train_samples), "samples\n")
  cat("Test set created with", length(test_samples), "samples\n")
  cat("Training subtypes:", table(train_data$subtypes$subtype), "\n")
  cat("Test subtypes:", table(test_data$subtypes$subtype), "\n")
  
  return(list(train = train_data, test = test_data))
}

# Main execution
cat("=== Multi-omics Data Download ===\n")

# Download all data types
mrna_data <- download_mrna_data()
mirna_data <- download_mirna_data()
protein_data <- download_protein_data()
subtype_info <- get_subtype_info()

# Create train/test splits
datasets <- create_train_test_sets(mrna_data, mirna_data, protein_data, subtype_info)

# Save summary information
download_summary <- list(
  mrna_samples = ncol(mrna_data$counts),
  mirna_samples = ncol(mirna_data$counts),
  protein_samples = ncol(protein_data$matrix),
  train_samples = ncol(datasets$train$mrna),
  test_samples = ncol(datasets$test$mrna),
  subtypes = table(subtype_info$subtype)
)

save(download_summary, file = file.path(multiomics_dir, "download_summary.RData"))

cat("Multi-omics data download completed successfully!\n")
cat("Data saved to:", multiomics_dir, "\n")