# Microarray Analysis Pipeline Documentation

## Overview
This repository contains the R code and documentation for microarray data analysis performed for [TGFβ/cyclin D1/Smad-mediated inhibition of BMP4 promotes breast cancer stem cell self-renewal activity]. The pipeline includes differential expression analysis, GO term enrichment, and visualization of gene expression patterns across multiple biological processes.

## Requirements

### R Version
- R >= 4.0.0

### Required R Packages
```R
- limma
- pheatmap
- ggplot2
- ggrepel
- tidyverse
- dplyr
- GEOquery
- beadarray
- illuminaHumanv1.db
- illuminaHumanv2.db
- illuminaHumanv3.db
- BeadArrayUseCases
- GOstats
- GenomicRanges
- Biostrings
```

## Data Structure
### Input Files
1. `SCP2_Microarray.txt`
   - Raw microarray expression data
   - Tab-delimited text file
   - Row names: gene identifiers
   - Columns: samples

2. `GO_Biological_Process_2018_table.txt`
   - GO term annotations
   - Tab-delimited text file
   - Contains gene sets for biological processes

3. `BMPs.txt`
   - BMP family gene expression data
   - Tab-delimited text file

## Analysis Pipeline

### 1. Data Preprocessing
- Log2 transformation of raw data
- Quantile normalization
- Quality control using boxplots

### 2. Differential Expression Analysis
- Linear model fitting using limma
- Empirical Bayes statistics
- Multiple testing correction using Benjamini-Hochberg method

### 3. GO Term Analysis
Analyzes gene expression patterns in multiple biological processes:
- Cell migration
- Cell proliferation
- Cell differentiation
- Cell adhesion
- Extracellular matrix organization
- Blood coagulation
- Signal transduction
- Endoderm development
- Peptidyl-tyrosine modification

### 4. Visualization
- Correlation heatmaps
- Volcano plots
- Process-specific heatmaps

## Output Files
1. `microarray.csv`: Differential expression results
2. `normalized_data.csv`: Normalized expression data
3. `volcano_plot.png/jpeg`: Volcano plot visualization
4. Multiple heatmap visualizations for each biological process

## Usage
1. Install required packages:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('beadarray', 'GEOquery', 'illuminaHumanv1.db',
                      'illuminaHumanv2.db', 'illuminaHumanv3.db', 
                      'BeadArrayUseCases', 'GOstats', 'GenomicRanges',
                      'Biostrings'))
```

2. Place input files in the working directory
3. Run the analysis script
4. Check output files for results

## Methods
### Normalization
Quantile normalization was performed on log2-transformed expression data to remove systematic variations and make samples comparable.

### Differential Expression
Empirical Bayes moderated t-statistics were used to identify differentially expressed genes. P-values were adjusted for multiple testing using the Benjamini-Hochberg method.

### GO Term Analysis
Gene sets were analyzed based on Gene Ontology biological process terms, focusing on specific cellular processes relevant to [specific biological context].

## Citation
If you use this pipeline, please cite:

[TGFβ/cyclin D1/Smad-mediated inhibition of BMP4 promotes breast cancer stem cell self-renewal activity
G. Yan, M. Dai, C. Zhang, S. Poulet, A. Moamer, N. Wang, et al.
Oncogenesis 2021 Vol. 10 Issue 3 Pages 21
DOI: 10.1038/s41389-021-00310-5
]

## Contact
[gang.yan@mail.mcgill.ca]

## Additional Notes
- Analysis parameters can be modified in the script as needed
- Default significance thresholds: adjusted p-value < 0.05
- Heatmaps are row-scaled for better visualization of expression patterns


