# Multi-omics Data Download Script
# Downloads TCGA breast cancer data (mRNA, miRNA, protein) following DIABLO methodology

suppressMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(SummarizedExperiment)
})

cat("Loading multi-omics data download functions...\n")

# Create multi-omics data directory
multiomics_dir <- file.path(data_dir, "multiomics")
dir.create(multiomics_dir, recursive = TRUE, showWarnings = FALSE)

# Download mRNA data for breast cancer
download_mrna_data <- function() {
  cat("Downloading mRNA data for TCGA-BRCA...\n")
  
  # Query mRNA data
  query_mrna <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download the data
  GDCdownload(query_mrna, method = "api", files.per.chunk = 10)
  
  # Prepare the data
  mrna_data <- GDCprepare(query_mrna)
  
  # Extract counts matrix
  mrna_counts <- assay(mrna_data)
  
  # Extract sample information
  mrna_samples <- as.data.frame(colData(mrna_data))
  
  # Save data
  save(mrna_counts, mrna_samples, file = file.path(multiomics_dir, "mrna_data.RData"))
  cat("mRNA data downloaded and saved.\n")
  
  return(list(counts = mrna_counts, samples = mrna_samples))
}

# Download miRNA data for breast cancer
download_mirna_data <- function() {
  cat("Downloading miRNA data for TCGA-BRCA...\n")
  
  # Query miRNA data
  query_mirna <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification"
  )
  
  # Download the data
  GDCdownload(query_mirna, method = "api", files.per.chunk = 10)
  
  # Prepare the data
  mirna_data <- GDCprepare(query_mirna)
  
  # Extract counts matrix
  mirna_counts <- assay(mirna_data)
  
  # Extract sample information
  mirna_samples <- as.data.frame(colData(mirna_data))
  
  # Save data
  save(mirna_counts, mirna_samples, file = file.path(multiomics_dir, "mirna_data.RData"))
  cat("miRNA data downloaded and saved.\n")
  
  return(list(counts = mirna_counts, samples = mirna_samples))
}

# Download protein data for breast cancer
download_protein_data <- function() {
  cat("Downloading protein data for TCGA-BRCA...\n")
  
  # Query protein data
  query_protein <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification"
  )
  
  # Download the data
  GDCdownload(query_protein, method = "api", files.per.chunk = 10)
  
  # Prepare the data
  protein_data <- GDCprepare(query_protein)
  
  # Extract protein matrix
  protein_matrix <- assay(protein_data)
  
  # Extract sample information
  protein_samples <- as.data.frame(colData(protein_data))
  
  # Save data
  save(protein_matrix, protein_samples, file = file.path(multiomics_dir, "protein_data.RData"))
  cat("Protein data downloaded and saved.\n")
  
  return(list(matrix = protein_matrix, samples = protein_samples))
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