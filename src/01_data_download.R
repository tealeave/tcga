# TCGA Data Download Script
# Downloads LUAD and LUSC RNA-seq data using TCGAbiolinks

suppressMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

cat("Loading TCGA data download functions...\n")

# Download LUAD data
cat("Downloading LUAD data...\n")
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_luad, method = "api", files.per.chunk = 10)
luad_data <- GDCprepare(query_luad)

# Download LUSC data  
cat("Downloading LUSC data...\n")
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_lusc, method = "api", files.per.chunk = 10)
lusc_data <- GDCprepare(query_lusc)

# Save raw data
cat("Saving raw data...\n")
save(luad_data, file = file.path(data_dir, "luad_raw.RData"))
save(lusc_data, file = file.path(data_dir, "lusc_raw.RData"))

# Extract counts matrices
luad_counts <- assay(luad_data, "unstranded")
lusc_counts <- assay(lusc_data, "unstranded")

# Extract sample information
luad_samples <- colData(luad_data)
lusc_samples <- colData(lusc_data)

# Save processed counts and sample info
save(luad_counts, luad_samples, file = file.path(data_dir, "luad_counts.RData"))
save(lusc_counts, lusc_samples, file = file.path(data_dir, "lusc_counts.RData"))

cat("Data download completed successfully!\n")
cat("LUAD samples:", ncol(luad_counts), "\n")
cat("LUSC samples:", ncol(lusc_counts), "\n")
cat("Genes:", nrow(luad_counts), "\n")