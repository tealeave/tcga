# Data Preprocessing Script
# Gene ID conversion and matrix combination following blog methodology

suppressMessages({
  library(org.Hs.eg.db)
  library(dplyr)
})

cat("Loading preprocessing functions...\n")

# Load downloaded data
load(file.path(data_dir, "luad_counts.RData"))
load(file.path(data_dir, "lusc_counts.RData"))

# Function to convert ENSEMBL IDs to gene symbols
convert_ensembl_to_symbols <- function(count_matrix) {
  cat("Converting ENSEMBL IDs to gene symbols...\n")
  
  # Extract ENSEMBL IDs (remove version numbers)
  ensembl_ids <- sub("\\..*", "", rownames(count_matrix))
  
  # Convert to gene symbols
  symbols <- mapIds(org.Hs.eg.db, 
                   keys = ensembl_ids,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")
  
  # Remove genes without symbols
  valid_symbols <- !is.na(symbols)
  count_matrix <- count_matrix[valid_symbols, ]
  symbols <- symbols[valid_symbols]
  
  # Remove duplicated gene symbols (keep first occurrence)
  unique_symbols <- !duplicated(symbols)
  count_matrix <- count_matrix[unique_symbols, ]
  symbols <- symbols[unique_symbols]
  
  # Set gene symbols as rownames
  rownames(count_matrix) <- symbols
  
  cat("Converted", length(symbols), "genes\n")
  return(count_matrix)
}

# Convert gene IDs for both datasets
luad_counts_symbols <- convert_ensembl_to_symbols(luad_counts)
lusc_counts_symbols <- convert_ensembl_to_symbols(lusc_counts)

# Find common genes between datasets
common_genes <- intersect(rownames(luad_counts_symbols), rownames(lusc_counts_symbols))
cat("Found", length(common_genes), "common genes\n")

# Subset to common genes
luad_counts_common <- luad_counts_symbols[common_genes, ]
lusc_counts_common <- lusc_counts_symbols[common_genes, ]

# Create sample type labels
luad_labels <- rep("LUAD", ncol(luad_counts_common))
lusc_labels <- rep("LUSC", ncol(lusc_counts_common))

# Combine matrices
combined_counts <- cbind(luad_counts_common, lusc_counts_common)
sample_labels <- c(luad_labels, lusc_labels)

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = colnames(combined_counts),
  cancer_type = sample_labels,
  stringsAsFactors = FALSE
)

# Save processed data
save(combined_counts, sample_metadata, file = file.path(data_dir, "combined_counts.RData"))

# Save individual processed datasets
save(luad_counts_symbols, lusc_counts_symbols, file = file.path(data_dir, "processed_counts.RData"))

cat("Data preprocessing completed successfully!\n")
cat("Combined matrix dimensions:", nrow(combined_counts), "x", ncol(combined_counts), "\n")
cat("LUAD samples:", sum(sample_labels == "LUAD"), "\n")
cat("LUSC samples:", sum(sample_labels == "LUSC"), "\n")