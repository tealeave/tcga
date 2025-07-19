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
  
  # Extract TCGA sample barcodes from cases field
  log_info("Extracting TCGA sample barcodes from cases field...")
  
  if ("cases" %in% names(files_info)) {
    # Get TCGA barcodes directly from cases
    tcga_barcodes <- sapply(files_info$cases, function(x) {
      if (length(x) > 0) {
        return(as.character(x[[1]]))
      }
      return(NA)
    })
    
    # Filter valid barcodes
    valid_idx <- !is.na(tcga_barcodes)
    tcga_barcodes <- tcga_barcodes[valid_idx]
    
    # Filter counts matrix to match valid barcodes
    mirna_counts <- mirna_counts[, valid_idx]
    
    # Use TCGA barcodes directly as sample identifiers
    colnames(mirna_counts) <- tcga_barcodes
    
    # Create sample metadata with TCGA barcodes
    mirna_samples <- data.frame(
      submitter_id = tcga_barcodes,
      sample_type = "Primary Tumor",
      stringsAsFactors = FALSE
    )
    rownames(mirna_samples) <- tcga_barcodes
    
  } else {
    # Fallback: use file IDs (should not happen with proper query)
    log_warn("Could not extract TCGA barcodes, using file IDs")
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
  
  # Find all files and directories matching the pattern in current directory
  all_items <- list.files(".", full.names = FALSE, include.dirs = TRUE)
  temp_items <- all_items[grepl(temp_pattern, all_items)]
  
  # Also check for .tar.gz files that might be left over
  temp_tar_files <- list.files(".", pattern = ".*_[0-9]+\\.tar\\.gz$", full.names = FALSE)
  
  all_temp_items <- c(temp_items, temp_tar_files)
  
  if (length(all_temp_items) > 0) {
    log_info("Found %d temporary items to clean up", length(all_temp_items))
    
    cleaned_count <- 0
    for (temp_item in all_temp_items) {
      if (file.exists(temp_item)) {
        # Check if it's a directory or file
        if (file.info(temp_item)$isdir) {
          unlink(temp_item, recursive = TRUE, force = TRUE)
          if (!file.exists(temp_item)) cleaned_count <- cleaned_count + 1
        } else {
          file.remove(temp_item)
          if (!file.exists(temp_item)) cleaned_count <- cleaned_count + 1
        }
      }
    }
    
    log_info("Successfully cleaned up %d temporary items", cleaned_count)
  } else {
    log_info("No temporary files found to clean up")
  }
  
  # Also clean up any .tmp files in the current directory
  tmp_files <- list.files(".", pattern = "\\.tmp$", full.names = FALSE)
  if (length(tmp_files) > 0) {
    file.remove(tmp_files)
    log_info("Cleaned up %d .tmp files", length(tmp_files))
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



# Fixed get_subtype_info function
get_subtype_info <- function() {
  operation_id <- log_operation_start("Subtype Data Download", "TCGA-BRCA Subtype Information")
  
  # Check if subtype data already exists
  save_path <- file.path(multiomics_dir, "subtype_info.RData")
  if (file.exists(save_path)) {
    log_info("Found existing subtype data, loading from file...")
    load(save_path)
    log_data_summary(subtype_info, "subtype information (loaded)", "data.frame")
    log_operation_end(operation_id, "Success (loaded existing data)")
    return(subtype_info)
  }

  result <- with_error_handling({
    # Use the dedicated function to get curated subtype data
    log_info("Downloading subtype data from PanCancerAtlas...")
    tcga_subtypes <- TCGAbiolinks::PanCancerAtlas_subtypes()
    
    # First, let's examine the structure to understand available columns
    log_info("Available columns in subtype data: %s", paste(colnames(tcga_subtypes), collapse = ", "))
    
    # Filter for breast cancer (BRCA) first
    brca_subtypes <- tcga_subtypes[tcga_subtypes$cancer.type == "BRCA", ]
    log_info("Found %d BRCA samples in subtype data", nrow(brca_subtypes))
    
    # Check available columns and select appropriate ones
    available_cols <- colnames(brca_subtypes)
    
    # Look for sample ID column (common names)
    sample_id_col <- NULL
    for (col_name in c("pan.samplesID", "sample", "sampleID", "submitter_id")) {
      if (col_name %in% available_cols) {
        sample_id_col <- col_name
        break
      }
    }
    
    # Look for subtype column (common names)
    subtype_col <- NULL
    for (col_name in c("Subtype_Selected", "subtype", "SUBTYPE", "molecular_subtype")) {
      if (col_name %in% available_cols) {
        subtype_col <- col_name
        break
      }
    }
    
    if (is.null(sample_id_col) || is.null(subtype_col)) {
      log_warn("Could not find expected columns. Available columns: %s", paste(available_cols, collapse = ", "))
      stop("Required columns not found in subtype data")
    }
    
    log_info("Using sample ID column: %s", sample_id_col)
    log_info("Using subtype column: %s", subtype_col)
    
    # Select and rename columns using standard R syntax
    subtype_info <- brca_subtypes[, c(sample_id_col, subtype_col)]
    colnames(subtype_info) <- c("patient_id", "subtype")
    
    # Filter out missing subtypes
    subtype_info <- subtype_info[!is.na(subtype_info$subtype) & subtype_info$subtype != "", ]
    
    # Convert patient IDs to 12-character format (TCGA-XX-XXXX)
    subtype_info$patient_id <- substr(subtype_info$patient_id, 1, 12)
    
    # Remove duplicate patient IDs (keep first occurrence)
    subtype_info <- subtype_info[!duplicated(subtype_info$patient_id), ]
    
    # Convert to data.frame and set row names
    subtype_info <- as.data.frame(subtype_info)
    rownames(subtype_info) <- subtype_info$patient_id
    
    log_data_summary(subtype_info, "subtype information", "data.frame")
    log_info("Subtype distribution: %s", paste(names(table(subtype_info$subtype)), "=", table(subtype_info$subtype), collapse = ", "))
    
    # Save the cleaned subtype information
    save(subtype_info, file = save_path)
    log_file_operation("saved", save_path, TRUE, "subtype data")
    
    subtype_info
  }, "Subtype Data Download")
  
  log_operation_end(operation_id, if(!is.null(result)) "Success" else "Failed")
  return(result)
}


# Create training and test sets following DIABLO methodology
create_train_test_sets <- function(mrna_data, mirna_data, protein_data, subtype_info) {
  cat("Creating training and test sets...\n")
  
  # Log data dimensions after loading
  log_info("Data loaded - mRNA: %dx%d, miRNA: %dx%d, protein: %dx%d", 
           nrow(mrna_data$counts), ncol(mrna_data$counts),
           nrow(mirna_data$counts), ncol(mirna_data$counts), 
           nrow(protein_data$matrix), ncol(protein_data$matrix))
  
  # Extract patient barcodes from sample barcodes (use exact patient ID extraction)
  
  # Log actual sample names to understand naming patterns
  log_info("Sample mRNA names (first 5): %s", paste(head(colnames(mrna_data$counts), 5), collapse=", "))
  log_info("Sample miRNA names (first 5): %s", paste(head(colnames(mirna_data$counts), 5), collapse=", "))
  log_info("Sample protein names (first 5): %s", paste(head(colnames(protein_data$matrix), 5), collapse=", "))
  
  # Use consistent patient ID extraction across all datasets
  extract_patient_id <- function(sample_barcode) {
    # Extract the patient ID (first 12 characters after TCGA-)
    if (grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_barcode)) {
      return(substr(sample_barcode, 1, 12))
    } else {
      return(sample_barcode)  # fallback
    }
  }
  
  mrna_patients <- sapply(colnames(mrna_data$counts), extract_patient_id)
  mirna_patients <- sapply(colnames(mirna_data$counts), extract_patient_id)
  protein_patients <- sapply(colnames(protein_data$matrix), extract_patient_id)
  
  # Log extracted patient IDs to verify extraction logic
  log_info("Extracted mRNA patient IDs (first 5): %s", paste(head(unique(mrna_patients), 5), collapse=", "))
  log_info("Extracted miRNA patient IDs (first 5): %s", paste(head(unique(mirna_patients), 5), collapse=", "))
  log_info("Extracted protein patient IDs (first 5): %s", paste(head(unique(protein_patients), 5), collapse=", "))
  
  # Get common patients across all data types
  common_patients <- intersect(intersect(mrna_patients, mirna_patients), protein_patients)
  
  # Log patient extraction results
  log_info("Patients extracted - mRNA: %d, miRNA: %d, protein: %d, common: %d", 
           length(unique(mrna_patients)), length(unique(mirna_patients)), 
           length(unique(protein_patients)), length(common_patients))

  # Find actual sample barcodes that exist across all three datasets for common patients
  # Create mapping from patient IDs to sample barcodes for each omics type
  mrna_patient_sample_map <- data.frame(
    patient_id = mrna_patients,
    sample_id = colnames(mrna_data$counts),
    stringsAsFactors = FALSE
  )

  # Print the first few rows of the mapping for verification
  log_info("mRNA patient-sample mapping (first 5 rows):\n%s", 
           paste(capture.output(head(mrna_patient_sample_map)), collapse="\n"))
  
  mirna_patient_sample_map <- data.frame(
    patient_id = mirna_patients,
    sample_id = colnames(mirna_data$counts),
    stringsAsFactors = FALSE
  )
  
  protein_patient_sample_map <- data.frame(
    patient_id = protein_patients,
    sample_id = colnames(protein_data$matrix),
    stringsAsFactors = FALSE
  )
  
  # Find patients who have samples in all three datasets
  patients_with_all_omics <- common_patients[
    common_patients %in% mrna_patient_sample_map$patient_id &
    common_patients %in% mirna_patient_sample_map$patient_id &
    common_patients %in% protein_patient_sample_map$patient_id
  ]
  
  # Create unified sample sets for each patient across all omics, this is a df with 1 col "patien_id"
  unified_samples <- data.frame(
    patient_id = patients_with_all_omics,
    stringsAsFactors = FALSE
  )
  
  # Map each patient to their actual sample barcodes in each omics
  unified_samples$mrna_sample <- sapply(patients_with_all_omics, function(patient) {
    mrna_patient_sample_map$sample_id[mrna_patient_sample_map$patient_id == patient][1]
  })
  
  unified_samples$mirna_sample <- sapply(patients_with_all_omics, function(patient) {
    mirna_patient_sample_map$sample_id[mirna_patient_sample_map$patient_id == patient][1]
  })
  
  unified_samples$protein_sample <- sapply(patients_with_all_omics, function(patient) {
    protein_patient_sample_map$sample_id[protein_patient_sample_map$patient_id == patient][1]
  })
  
  log_info("Unified samples created for %d patients across all omics", nrow(unified_samples))

  # Print the first few rows of unified samples for verification
  log_info("Unified samples (first 5 rows):\n%s", 
           paste(capture.output(head(unified_samples)), collapse="\n"))

  # Filter subtype info for patients with all omics data
  available_subtypes <- subtype_info[subtype_info$patient_id %in% unified_samples$patient_id, ]
  
  # Log available subtypes before filtering
  log_info("Available subtypes: %s", paste(unique(available_subtypes$subtype), collapse=", "))
  log_info("Patients with subtype info: %d", nrow(available_subtypes))
  
  # Focus on main subtypes: Basal, Her2, LumA, LumB (exclude Normal)
  main_subtypes <- available_subtypes[!available_subtypes$subtype %in% c("BRCA.Normal"), ]
  
  # Log filtering results
  log_info("Patients after main subtype filtering: %d", nrow(main_subtypes))
  
  # Create training set (150 patients, 50 per subtype)
  set.seed(42)  # For reproducibility
  train_patients <- main_subtypes %>%
    group_by(subtype) %>%
    sample_n(min(50, n())) %>%  # 50 patients per subtype, or all if fewer
    ungroup() %>%
    pull(patient_id)
  
  # Filter unified_samples for training patients and use omics-specific samples for subsetting
  train_unified <- unified_samples[unified_samples$patient_id %in% train_patients, ]
  
  # Create training data using omics-specific sample barcodes
  train_data <- list(
    mrna = mrna_data$counts[, train_unified$mrna_sample],
    mirna = mirna_data$counts[, train_unified$mirna_sample],
    protein = protein_data$matrix[, train_unified$protein_sample],
    subtypes = data.frame(
      sample_id = train_unified$mrna_sample,  # Use mRNA samples as reference for subtypes
      subtype = main_subtypes$subtype[match(train_unified$patient_id, main_subtypes$patient_id)]
    )
  )
  
  # Create test set (up to 75 patients, ~25 per subtype, with only mRNA and miRNA)
  remaining_patients <- setdiff(main_subtypes$patient_id, train_patients)
  test_patients <- main_subtypes %>%
    filter(patient_id %in% remaining_patients) %>%
    group_by(subtype) %>%
    sample_n(min(25, n())) %>%
    ungroup() %>%
    pull(patient_id)
  
  # Filter unified_samples for test patients
  test_unified <- unified_samples[unified_samples$patient_id %in% test_patients, ]
  
  # Create test data using omics-specific sample barcodes (protein missing as per methodology)
  test_data <- list(
    mrna = mrna_data$counts[, test_unified$mrna_sample],
    mirna = mirna_data$counts[, test_unified$mirna_sample],
    protein = NULL,  # Missing protein data as in the DIABLO methodology
    subtypes = data.frame(
      sample_id = test_unified$mrna_sample,  # Use mRNA samples as reference for subtypes
      subtype = main_subtypes$subtype[match(test_unified$patient_id, main_subtypes$patient_id)]
    )
  )
  
  # Save training and test sets
  save(train_data, file = file.path(multiomics_dir, "train_data.RData"))
  save(test_data, file = file.path(multiomics_dir, "test_data.RData"))
  
  cat("Training set created with", nrow(train_unified), "samples (", length(train_patients), "patients)\n")
  cat("Test set created with", nrow(test_unified), "samples (", length(test_patients), "patients)\n")
  
  # Log final results
  log_info("Created training set: %d samples (%d patients)", nrow(train_unified), length(train_patients))
  log_info("Created test set: %d samples (%d patients)", nrow(test_unified), length(test_patients))
  cat("Training subtypes:", table(train_data$subtypes$subtype), "\n")
  cat("Test subtypes:", table(test_data$subtypes$subtype), "\n")
  
  return(list(train = train_data, test = test_data))
}


# Main execution
# Ensure logs directory exists and setup output capture
if (!dir.exists("logs")) {
  dir.create("logs", showWarnings = FALSE)
}

# Simple global output capture - only capture if not already open
tryCatch({
  # Check if sink is already active
  if (sink.number() == 0) {
    sink("logs/tcga_pipeline_output.log", append = TRUE)
    sink("logs/tcga_pipeline_output.log", type = "message", append = TRUE)
  }
}, error = function(e) {
  warning("Could not set up output capture: ", e$message)
})

cat("=== Multi-omics Data Download ===\n")

# Download all data types
mrna_data <- download_mrna_data()
mirna_data <- download_mirna_data()
protein_data <- download_protein_data()
subtype_info <- get_subtype_info()

# Print a couple rows of the subtype information for verification
log_info("Subtype information (first 5 rows):\n%s",
         paste(capture.output(head(subtype_info)), collapse="\n"))

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

# DEBUG: Print detailed data structure information
log_info("=== DEBUG: Final Data Structure Analysis ===")

# Log training data structure
log_info("--- Training Data Heads ---")
for (block_name in c("mrna", "mirna", "protein")) {
  if (!is.null(datasets$train[[block_name]])) {
    log_data_summary(datasets$train[[block_name]], paste(block_name, "training data"), "matrix")
    log_info("First 5 rows x 5 cols:\n%s", 
             paste(capture.output(head(datasets$train[[block_name]][1:5, 1:5])), collapse="\n"))
    log_info("Column names (samples, first 10): %s", 
             paste(head(colnames(datasets$train[[block_name]]), 10), collapse=", "))
  }
}

log_data_summary(datasets$train$subtypes, "Training subtypes", "data.frame")

# Log test data structure  
log_info("--- Test Data Heads ---")
for (block_name in c("mrna", "mirna", "protein")) {
  if (!is.null(datasets$test[[block_name]])) {
    log_data_summary(datasets$test[[block_name]], paste(block_name, "test data"), "matrix")
    log_info("First 5 rows x 5 cols:\n%s", 
             paste(capture.output(head(datasets$test[[block_name]][1:5, 1:5])), collapse="\n"))
    log_info("Column names (samples, first 10): %s", 
             paste(head(colnames(datasets$test[[block_name]]), 10), collapse=", "))
  }
}

log_data_summary(datasets$test$subtypes, "Test subtypes", "data.frame")

# Log sample counts
log_info("Training sample counts - mRNA: %d, miRNA: %d, protein: %d, subtypes: %d", 
         ncol(datasets$train$mrna), ncol(datasets$train$mirna), 
         ncol(datasets$train$protein), nrow(datasets$train$subtypes))

log_info("Test sample counts - mRNA: %d, miRNA: %d, subtypes: %d", 
         ncol(datasets$test$mrna), ncol(datasets$test$mirna), nrow(datasets$test$subtypes))

# Check sample name matching
train_mrna_match <- all(colnames(datasets$train$mrna) == datasets$train$subtypes$sample_id)
train_mirna_match <- all(colnames(datasets$train$mirna) == datasets$train$subtypes$sample_id)
train_protein_match <- all(colnames(datasets$train$protein) == datasets$train$subtypes$sample_id)

test_mrna_match <- all(colnames(datasets$test$mrna) == datasets$test$subtypes$sample_id)
test_mirna_match <- all(colnames(datasets$test$mirna) == datasets$test$subtypes$sample_id)

log_info("Training set name matching - mRNA vs subtypes: %s, miRNA vs subtypes: %s, protein vs subtypes: %s", 
         train_mrna_match, train_mirna_match, train_protein_match)

log_info("Test set name matching - mRNA vs subtypes: %s, miRNA vs subtypes: %s", 
         test_mrna_match, test_mirna_match)

cat("Multi-omics data download completed successfully!\n")
cat("Data saved to:", multiomics_dir, "\n")

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