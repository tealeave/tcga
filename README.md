# TCGA Multi-omics Analysis Pipeline

This project implements comprehensive analysis pipelines for TCGA data, including both traditional PCA analysis and advanced multi-omics integration using DIABLO (Data Integration Analysis for Biomarker discovery using Latent variable approaches).

## Analysis Types

### 1. PCA Analysis (Original)
- LUAD (Lung Adenocarcinoma) and LUSC (Lung Squamous Cell Carcinoma) samples
- Principal Component Analysis (PCA) methodology
- Follows "Diving into Genetics and Genomics" blog posts

### 2. Multi-omics DIABLO Analysis (New)
- TCGA breast cancer data (BRCA)
- Multi-omics integration: mRNA, miRNA, and protein data
- Breast cancer subtype classification (Basal, Her2, LumA)
- Follows mixOmics DIABLO methodology

## Project Structure

```
tcga/
├── src/                    # R analysis scripts
│   ├── 01_data_download.R      # TCGA data download using TCGAbiolinks
│   ├── 02_data_preprocessing.R # Gene ID conversion and matrix combination
│   ├── 03_pca_analysis.R       # PCA analysis (raw counts vs CPM)
│   ├── 04_visualization.R      # Generate all plots from blog posts
│   ├── 05_marker_analysis.R    # Marker gene analysis
│   ├── 06_multiomics_download.R    # Multi-omics data download
│   ├── 07_multiomics_preprocessing.R # Multi-omics data preprocessing
│   ├── 08_block_splsda.R       # Block sPLS-DA analysis
│   ├── 09_diablo_analysis.R    # DIABLO analysis pipeline
│   └── 10_multiomics_visualization.R # Multi-omics visualization suite
├── data/                   # Raw and processed data (gitignored)
│   └── multiomics/        # Multi-omics specific data
├── results/               # Analysis outputs (gitignored)
│   ├── plots/             # Generated visualizations
│   └── tables/            # Analysis result tables
├── docs/                  # Documentation
├── pipeline.R             # Main driver script (PCA analysis)
├── pipeline_multiomics.R  # Multi-omics pipeline driver
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

### PCA Analysis (Original)

To run the complete PCA analysis pipeline:

```bash
# Load R module
module load R/4.4.2

# Run the full PCA pipeline
Rscript pipeline.R
```

### Multi-omics DIABLO Analysis (New)

To run the multi-omics DIABLO analysis:

```bash
# Load R module
module load R/4.4.2

# Run the multi-omics pipeline
Rscript pipeline_multiomics.R
```

### Individual Steps

You can also run individual analysis steps:

**PCA Analysis:**
```r
# In R console
source("src/01_data_download.R")      # Download TCGA data
source("src/02_data_preprocessing.R") # Process and combine data
source("src/03_pca_analysis.R")       # Perform PCA analysis
source("src/04_visualization.R")      # Generate plots
source("src/05_marker_analysis.R")    # Analyze marker genes
```

**Multi-omics Analysis:**
```r
# In R console
source("src/06_multiomics_download.R")    # Download multi-omics data
source("src/07_multiomics_preprocessing.R") # Preprocess multi-omics data
source("src/08_block_splsda.R")           # Block sPLS-DA analysis
source("src/09_diablo_analysis.R")        # DIABLO analysis
source("src/10_multiomics_visualization.R") # Generate visualizations
```

## Analysis Overview

### PCA Analysis (Original Pipeline)

#### 1. Data Download (`01_data_download.R`)
- Downloads LUAD and LUSC RNA-seq data from TCGA using TCGAbiolinks
- Extracts raw counts matrices and sample metadata
- Saves processed data for downstream analysis

#### 2. Data Preprocessing (`02_data_preprocessing.R`)
- Converts ENSEMBL gene IDs to gene symbols using org.Hs.eg.db
- Combines LUAD and LUSC datasets
- Filters for common genes between datasets
- Creates sample metadata with cancer type labels

#### 3. PCA Analysis (`03_pca_analysis.R`)
- Selects top 1000 most variable genes
- Performs PCA on:
  - Raw counts (shows sequencing depth correlation)
  - CPM-normalized counts (better cancer type separation)
- Calculates variance explained by each principal component

#### 4. Visualization (`04_visualization.R`)
- Generates PCA scatter plots for both raw and normalized data
- Creates variance explained plots
- Produces heatmaps of top variable genes
- Saves all plots in multiple formats (PNG, PDF)

#### 5. Marker Gene Analysis (`05_marker_analysis.R`)
- Analyzes expression of lung cancer marker genes:
  - NAPSA (lung adenocarcinoma marker)
  - TFF1 (lung adenocarcinoma marker)
  - TP63 (lung squamous cell carcinoma marker)
  - KRT5 (lung squamous cell carcinoma marker)
- Performs statistical testing between cancer types
- Creates expression distribution plots

### Multi-omics DIABLO Analysis (New Pipeline)

#### 1. Multi-omics Data Download (`06_multiomics_download.R`)
- Downloads TCGA breast cancer (BRCA) multi-omics data:
  - mRNA expression data
  - miRNA expression data
  - Protein expression data
  - Clinical subtype information
- Creates training set (150 samples) and test set (70 samples)
- Handles missing protein data in test set as per DIABLO methodology

#### 2. Multi-omics Preprocessing (`07_multiomics_preprocessing.R`)
- Standardizes data matrices across omics types
- Performs feature selection:
  - mRNA: Top 5000 most variable genes
  - miRNA: Top 200 most variable miRNAs
  - Protein: Filters and imputes missing values
- Handles sample matching across data types
- Creates design matrix for DIABLO analysis

#### 3. Block sPLS-DA Analysis (`08_block_splsda.R`)
- Implements multiblock sparse Partial Least Squares Discriminant Analysis
- Tunes optimal number of components via cross-validation
- Optimizes variable selection (keepX parameters)
- Builds final multiblock model
- Evaluates model performance

#### 4. DIABLO Analysis (`09_diablo_analysis.R`)
- Implements full DIABLO workflow
- Handles missing data blocks (protein data in test set)
- Generates predictions using multiple methods:
  - Centroids distance
  - Majority vote
  - Weighted vote
- Extracts biomarker signatures
- Assesses model stability through bootstrap analysis

#### 5. Multi-omics Visualization (`10_multiomics_visualization.R`)
- Creates comprehensive visualization suite:
  - Sample plots (component scatter plots)
  - Variable correlation networks
  - Clustered image maps (heatmaps)
  - Performance comparison plots
  - Biomarker importance plots
  - Arrow plots showing component contributions
  - Block correlation heatmaps
- Generates visualizations for both DIABLO and Block sPLS-DA results

## Key Findings

### PCA Analysis Findings

Based on the blog post methodology, the PCA analysis reveals:

1. **Raw counts PCA**: First principal component strongly correlates with sequencing depth
2. **CPM-normalized PCA**: Better separation between LUAD and LUSC samples
3. **Sample heterogeneity**: Some samples show unexpected expression patterns
4. **Marker genes**: Clear differential expression between cancer types

### Multi-omics DIABLO Analysis Findings

Based on the mixOmics DIABLO methodology, the multi-omics analysis reveals:

1. **Multi-omics integration**: Successful integration of mRNA, miRNA, and protein data
2. **Subtype classification**: Accurate classification of breast cancer subtypes (Basal, Her2, LumA)
3. **Biomarker discovery**: Identification of key biomarkers across different omics types
4. **Model performance**: High accuracy in cross-validation and test set predictions
5. **Missing data handling**: Robust performance even with missing protein data in test set
6. **Feature correlations**: Strong correlations between features across different omics blocks
7. **Stability assessment**: Consistent biomarker selection across bootstrap iterations

## Generated Outputs

### PCA Analysis Plots (`results/plots/`)
- `raw_counts_pca.png` - PCA plot using raw counts
- `cpm_normalized_pca.png` - PCA plot using CPM-normalized data
- `pca_comparison.png` - Side-by-side comparison
- `variance_explained.png` - Variance explained by each PC
- `top_genes_heatmap.png/pdf` - Heatmap of top variable genes
- `marker_genes_boxplot.png` - Marker gene expression by cancer type
- `marker_genes_violin.png` - Distribution plots for marker genes
- `pc1_scores_by_cancer_type.png` - PC1 scores visualization

### Multi-omics DIABLO Plots (`results/plots/`)
- `diablo_sample_plot_comp*.png` - Sample plots for different component combinations
- `diablo_network_comp*.png` - Variable correlation networks for each component
- `diablo_heatmap_comp*.png` - Clustered image maps for each component
- `diablo_component_tuning.png` - Component selection optimization
- `diablo_performance_comparison.png` - Prediction method comparison
- `diablo_biomarkers_*.png` - Biomarker importance plots for each data type
- `diablo_arrow_plot_comp*.png` - Arrow plots showing component contributions
- `diablo_correlation_*.png` - Block correlation heatmaps
- `diablo_stability_*.png` - Feature stability plots
- `diablo_model_configuration.png` - Model configuration summary
- `diablo_explained_variance.png` - Explained variance by component
- `block_spls-da_*.png` - Block sPLS-DA specific visualizations

### PCA Analysis Tables (`results/tables/`)
- `raw_pca_scores.csv` - PCA scores for raw counts analysis
- `cpm_pca_scores.csv` - PCA scores for CPM-normalized analysis
- `variance_explained.csv` - Variance explained by principal components
- `marker_gene_expression.csv` - Marker gene expression data
- `marker_gene_summary.csv` - Summary statistics for marker genes
- `marker_gene_tests.csv` - Statistical test results

### Multi-omics DIABLO Tables (`results/tables/`)
- `diablo_performance_summary.csv` - Model performance metrics
- `diablo_configuration.csv` - Model configuration parameters
- `diablo_biomarkers_*.csv` - Biomarker lists for each data type
- `diablo_confusion_*.csv` - Confusion matrices for different prediction methods
- `block_splsda_model_summary.csv` - Block sPLS-DA model summary
- `block_splsda_performance.csv` - Block sPLS-DA performance metrics
- `block_splsda_important_vars_*.csv` - Important variables for each data type

## Configuration

Analysis parameters can be modified in `config.yaml`:

### PCA Analysis Parameters
- Number of top variable genes for PCA
- Marker genes to analyze
- Plot dimensions and colors
- Statistical analysis parameters

### Multi-omics Analysis Parameters
- DIABLO analysis settings (components, cross-validation, keepX values)
- Block sPLS-DA parameters
- Feature selection thresholds
- Multi-omics visualization settings
- Breast cancer subtype definitions
- Training/test split configuration

## References

### PCA Analysis References
- Blog posts: [Diving into Genetics and Genomics](https://divingintogeneticsandgenomics.com/)
- [PCA TCGA Analysis Part 1](https://divingintogeneticsandgenomics.com/post/pca-tcga/)
- [PCA TCGA Analysis Part 2](https://divingintogeneticsandgenomics.com/post/pca-tcga2/)

### Multi-omics DIABLO Analysis References
- [mixOmics DIABLO Vignette](https://mixomicsteam.github.io/mixOmics-Vignette/id_06.html#id_06:diablo)
- [mixOmics Documentation](http://mixomics.org/)
- Singh, A., Shannon, C.P., Gautier, B., Rohart, F., Vacher, M., Tebbutt, S.J. and Lê Cao, K.A. (2019). DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics*, 35(17), 3055-3062.
- Rohart, F., Gautier, B., Singh, A., & Lê Cao, K. A. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. *PLoS computational biology*, 13(11), e1005752.

## Dependencies

### R Packages (PCA Analysis)
- **TCGAbiolinks**: TCGA data download
- **org.Hs.eg.db**: Gene ID conversion
- **edgeR**: RNA-seq normalization
- **ggplot2**: Visualization
- **ComplexHeatmap**: Heatmap generation
- **ggfortify**: PCA visualization
- **dplyr, tidyr**: Data manipulation

### R Packages (Multi-omics Analysis)
- **mixOmics**: Core multi-omics integration and DIABLO analysis
- **igraph**: Network analysis and visualization
- **corrplot**: Correlation matrix visualization
- **circlize**: Circular plots
- **SummarizedExperiment**: Bioconductor data structures
- **gridExtra**: Advanced plot arrangements
- **RColorBrewer**: Color palettes
- **yaml**: Configuration file handling

### Python Packages (minimal)
- **jupyter**: For notebook exploration

## Notes

- Raw TCGA data files are large and excluded from version control
- Analysis runtime depends on data download speed and system resources
- All plots are saved in high resolution (300 DPI) for publication quality
- Statistical tests use Benjamini-Hochberg correction for multiple comparisons