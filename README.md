# Code-for-Single-Cell-Analysis
 
# Single-Cell Analysis of Type 2 Diabetes Islet Data

This repository contains code and analysis for single-cell RNA sequencing (scRNA-seq) data from islet cells in individuals with Type 2 diabetes (T2D). The analysis is performed using data from the PANDB database, focusing on clustering, pseudotime analysis, and visualization. This is an on going project. 

## Features

- **Raw Data download**: Islet cell data download from PANCB database.
- **Quality control and mapping** : quality control and mapping performed by fastp and salmon.
- **Count matrix Generation**: Performed using tx2gene in R environment.
- **Data Processing**: Preprocessing and normalization of scRNA-seq data.
- **Clustering**: Identification of cell clusters using t-SNE and Seurat.
- **Pseudotime Analysis**: Inference of cell trajectories and pseudotime using Slingshot.
- **Visualization**: Generation of t-SNE plots and pseudotime trajectory plots.

## Installation

To run the analysis, you need to have R and Linux/UNIX enviironment, the following pacakges are needed:

fastp, salmon, tximport, GenomicFeatures, bioconductor, rhdf5, DESeq2, SIngleCellExperiment, biomaRt, Seurat, cowplot,
dplyr, ggplot2, slingshot, monocle3, SeuratWrappers

## Results
The project generates various visualizations, including t-SNE plots and pseudotime trajectory plots, which can help in understanding the cellular heterogeneity and lineage relationships in islet cells of T2D patients.

## Exmaple Results
- Violin plots of genes of interest 
- PCA plots
- PCA illustrating key features that contribute to each principal component
- t-SNE cluster plots for cell identification
- t-SNE clusters highlighting genes of interest

## Script definition

- **pancdb_Final_26_T2D_scRNA**: Bash script mass downloading raw scRNA-seq data from ftp links
- **Cleaning+Alignment**: Cleaning and aligning all raw data
- **CleaningAlignentReference**: Main script used for downloading and processing large scale scRNA-seq data
- **Creating-tx2gene**: Generating count matrix and saving data for further processing
- **Exploring-Count-Matrix**: Main script for generating figures and understanding the heterogeneity of beta cell failure in type 2 diabetes

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your changes. Ensure your code adheres to the project's coding standards and includes appropriate documentation.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For questions or feedback, please contact jonanzule@gmail.com.








