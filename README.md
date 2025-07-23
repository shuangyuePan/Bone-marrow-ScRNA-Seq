# Bone Marrow Single-Cell RNA Sequencing Analysis

[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v4.0+-orange)](https://satijalab.org/seurat/)

## Overview

This repository contains the computational analysis code for single-cell RNA sequencing (scRNA-seq) data from bone marrow samples, investigating the role of ARID1A in myeloid cell differentiation and neutrophil development. The study compares ARID1A knockout (KO) versus normal control (NC) conditions to understand chromatin remodeling effects on hematopoiesis.

## Study Design

- **Experimental Groups**: ARID1A knockout (KO) vs Normal control (NC)
- **Cell Types Analyzed**: Early Myeloid, Neutrophils, pre-MDSC, PMN-MDSC
- **Key Genes**: Arg2, Sat1, Ly6g, Csf1r, Il1b
- **Analysis Methods**: Trajectory analysis, differential expression, ChIP-seq integration

## Repository Structure

```
Bone-marrow-ScRNA-Seq/
├── README.md                    # This file
├── monecle/
│   └── all.R                   # Monocle2 trajectory analysis (Fig 4F, 4G, S9A, S9E)
├── Fig3F/
│   └── Arg2_peak_plot.R        # ChIP-seq peak visualization for Arg2 promoter
├── Fig4D/
│   └── density_plot.R          # Cell density analysis for neutrophil populations
├── Fig4E/
│   └── Feautre_plot.R          # Feature plots for key marker genes
├── Fig4H/
│   └── CellRatio.R             # Cell type proportion analysis
├── Fig4I/
│   └── density_plot.R          # Density analysis for pre-MDSC populations
└── Fig4J/
    └── density_plot.R          # Density analysis for PMN-MDSC populations
```

## Key Findings

### 1. Trajectory Analysis (Monocle2)
- **Pseudotime reconstruction** of myeloid differentiation from early myeloid cells to mature neutrophils
- **Branch point analysis** revealing divergent differentiation paths
- **State-specific gene expression** patterns across 5 distinct cellular states
- **Differential trajectory dynamics** between KO and NC conditions

### 2. ChIP-seq Integration
- **Arg2 promoter analysis** showing 10-fold enrichment in ARID1A KO vs NC
- **Chromatin accessibility changes** correlating with gene expression alterations

### 3. Cell Population Dynamics
- **Density-based clustering** analysis of neutrophil, pre-MDSC, and PMN-MDSC populations
- **Statistical significance testing** using Fisher's exact test
- **Neighbor analysis** with k-nearest neighbor algorithms

### 4. Marker Gene Expression
- **Ly6g**: Neutrophil maturation marker
- **Csf1r**: Myeloid lineage commitment
- **Arg2**: Metabolic reprogramming in MDSCs
- **Il1b**: Inflammatory response

## Dependencies

### R Packages Required
```r
# Core single-cell analysis
library(Seurat)          # v4.0+
library(monocle)         # v2.0+

# Data manipulation and visualization
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(pheatmap)

# Specialized analysis
library(harmony)         # Batch correction
library(FNN)            # K-nearest neighbors
library(rtracklayer)    # ChIP-seq data handling
library(GenomicRanges)  # Genomic intervals

# Utilities
library(qs)             # Fast serialization
library(Matrix)
library(parallel)
```

## Usage

### Prerequisites
1. **R version ≥ 4.0.0**
2. **Required R packages** (see Dependencies section)
3. **Input data**: `bm.qs` file containing processed Seurat object
4. **ChIP-seq data**: BigWig files for ARID1A ChIP-seq (Fig3F only)

### Running the Analysis

1. **Trajectory Analysis**:
```bash
cd monecle/
Rscript all.R
```

2. **Individual Figure Generation**:
```bash
# For each figure panel
cd Fig4D/
Rscript density_plot.R

cd Fig4E/
Rscript Feautre_plot.R

# ... and so on for other figures
```

### Expected Outputs

- **PDF files**: High-resolution publication-ready figures
- **Trajectory plots**: Cell state transitions and pseudotime
- **Density plots**: Population distribution analysis
- **Feature plots**: Gene expression visualization
- **Statistical results**: Printed to console with significance testing

## Computational Methods

### Single-Cell Analysis Pipeline
1. **Quality Control**: Cell and gene filtering
2. **Normalization**: Log-normalization and scaling
3. **Dimensionality Reduction**: PCA, t-SNE
4. **Clustering**: Graph-based clustering
5. **Trajectory Inference**: Monocle2 DDRTree algorithm

### Statistical Analysis
- **Differential Expression**: Negative binomial testing
- **Trajectory Significance**: q-value < 0.001 for cell type, q-value < 0.01 for treatment
- **Population Comparison**: Fisher's exact test with simulation (B=10,000)
- **Neighbor Analysis**: Radius-based search (r=2.66) with k-NN validation

## Citation

If you use this code in your research, please cite:

```bibtex
@article{author2024arid1a,
  title={ARID1A regulates myeloid cell differentiation and neutrophil development in bone marrow},
  author={[Author Names]},
  journal={[Journal Name]},
  year={2024},
  doi={10.xxxx/xxxxxx}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the analysis or code:
- **Corresponding Author**: [Email]
- **Computational Analysis**: [Email]
- **Issues**: Please use the GitHub issue tracker

## Acknowledgments

- Seurat development team for single-cell analysis framework
- Monocle development team for trajectory inference methods
- R/Bioconductor community for genomics tools

---

**Note**: This repository contains the computational analysis code Only. Raw sequencing data and processed count matrices are available through [Data Repository/GEO Accession].