# TCGA RNA-seq PCA Analysis

This project implements a comprehensive analysis pipeline for TCGA RNA-seq data, specifically analyzing LUAD (Lung Adenocarcinoma) and LUSC (Lung Squamous Cell Carcinoma) samples using Principal Component Analysis (PCA). The analysis follows the methodology described in the blog posts from "Diving into Genetics and Genomics".

## Project Structure

```
tcga/
├── src/                    # R analysis scripts
│   ├── 01_data_download.R      # TCGA data download using TCGAbiolinks
│   ├── 02_data_preprocessing.R # Gene ID conversion and matrix combination
│   ├── 03_pca_analysis.R       # PCA analysis (raw counts vs CPM)
│   ├── 04_visualization.R      # Generate all plots from blog posts
│   └── 05_marker_analysis.R    # Marker gene analysis
├── data/                   # Raw and processed data (gitignored)
├── results/               # Analysis outputs (gitignored)
│   ├── plots/             # Generated visualizations
│   └── tables/            # Analysis result tables
├── docs/                  # Documentation
├── pipeline.R             # Main driver script
├── config.yaml            # Analysis configuration
├── pyproject.toml         # Python dependencies (minimal)
├── renv.lock              # R package dependencies
└── README.md              # This file
```

## Setup Instructions

### Prerequisites

- Access to an HPC system with R 4.4.2 or higher
- Internet connection for downloading TCGA data

### 1. Environment Setup

**Load R module (HPC):**
```bash
module load R/4.4.2
```

**Activate R environment:**
```r
# In R console
renv::activate()
renv::restore()
```

### 2. Python Environment (optional)

The project includes minimal Python dependencies for potential notebook exploration:

```bash
# Python environment with uv
uv sync
```

## Running the Analysis

### Full Pipeline

To run the complete analysis pipeline:

```bash
# Load R module
module load R/4.4.2

# Run the full pipeline
Rscript pipeline.R
```

### Individual Steps

You can also run individual analysis steps:

```r
# In R console
source("src/01_data_download.R")      # Download TCGA data
source("src/02_data_preprocessing.R") # Process and combine data
source("src/03_pca_analysis.R")       # Perform PCA analysis
source("src/04_visualization.R")      # Generate plots
source("src/05_marker_analysis.R")    # Analyze marker genes
```

## Analysis Overview

### 1. Data Download (`01_data_download.R`)
- Downloads LUAD and LUSC RNA-seq data from TCGA using TCGAbiolinks
- Extracts raw counts matrices and sample metadata
- Saves processed data for downstream analysis

### 2. Data Preprocessing (`02_data_preprocessing.R`)
- Converts ENSEMBL gene IDs to gene symbols using org.Hs.eg.db
- Combines LUAD and LUSC datasets
- Filters for common genes between datasets
- Creates sample metadata with cancer type labels

### 3. PCA Analysis (`03_pca_analysis.R`)
- Selects top 1000 most variable genes
- Performs PCA on:
  - Raw counts (shows sequencing depth correlation)
  - CPM-normalized counts (better cancer type separation)
- Calculates variance explained by each principal component

### 4. Visualization (`04_visualization.R`)
- Generates PCA scatter plots for both raw and normalized data
- Creates variance explained plots
- Produces heatmaps of top variable genes
- Saves all plots in multiple formats (PNG, PDF)

### 5. Marker Gene Analysis (`05_marker_analysis.R`)
- Analyzes expression of lung cancer marker genes:
  - NAPSA (lung adenocarcinoma marker)
  - TFF1 (lung adenocarcinoma marker)
  - TP63 (lung squamous cell carcinoma marker)
  - KRT5 (lung squamous cell carcinoma marker)
- Performs statistical testing between cancer types
- Creates expression distribution plots

## Key Findings

Based on the blog post methodology, the analysis reveals:

1. **Raw counts PCA**: First principal component strongly correlates with sequencing depth
2. **CPM-normalized PCA**: Better separation between LUAD and LUSC samples
3. **Sample heterogeneity**: Some samples show unexpected expression patterns
4. **Marker genes**: Clear differential expression between cancer types

## Generated Outputs

### Plots (`results/plots/`)
- `raw_counts_pca.png` - PCA plot using raw counts
- `cpm_normalized_pca.png` - PCA plot using CPM-normalized data
- `pca_comparison.png` - Side-by-side comparison
- `variance_explained.png` - Variance explained by each PC
- `top_genes_heatmap.png/pdf` - Heatmap of top variable genes
- `marker_genes_boxplot.png` - Marker gene expression by cancer type
- `marker_genes_violin.png` - Distribution plots for marker genes
- `pc1_scores_by_cancer_type.png` - PC1 scores visualization

### Tables (`results/tables/`)
- `raw_pca_scores.csv` - PCA scores for raw counts analysis
- `cpm_pca_scores.csv` - PCA scores for CPM-normalized analysis
- `variance_explained.csv` - Variance explained by principal components
- `marker_gene_expression.csv` - Marker gene expression data
- `marker_gene_summary.csv` - Summary statistics for marker genes
- `marker_gene_tests.csv` - Statistical test results

## Configuration

Analysis parameters can be modified in `config.yaml`:

- Number of top variable genes for PCA
- Marker genes to analyze
- Plot dimensions and colors
- Statistical analysis parameters

## References

- Blog posts: [Diving into Genetics and Genomics](https://divingintogeneticsandgenomics.com/)
- [PCA TCGA Analysis Part 1](https://divingintogeneticsandgenomics.com/post/pca-tcga/)
- [PCA TCGA Analysis Part 2](https://divingintogeneticsandgenomics.com/post/pca-tcga2/)

## Dependencies

### R Packages
- **TCGAbiolinks**: TCGA data download
- **org.Hs.eg.db**: Gene ID conversion
- **edgeR**: RNA-seq normalization
- **ggplot2**: Visualization
- **ComplexHeatmap**: Heatmap generation
- **ggfortify**: PCA visualization
- **dplyr, tidyr**: Data manipulation

### Python Packages (minimal)
- **jupyter**: For notebook exploration

## Notes

- Raw TCGA data files are large and excluded from version control
- Analysis runtime depends on data download speed and system resources
- All plots are saved in high resolution (300 DPI) for publication quality
- Statistical tests use Benjamini-Hochberg correction for multiple comparisons