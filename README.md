# Code-for-Single-Cell-Analysis
 
# Single-Cell Analysis of Type 2 Diabetes Islet Data

This repository contains code and analysis for single-cell RNA sequencing (scRNA-seq) data from islet cells in individuals with Type 2 diabetes (T2D). The analysis is performed using data from the PANDb database, focusing on clustering, pseudotime analysis, and visualization.

## Features

- **Raw Data download**: Islet cell data download from PANCB database.
- **Quality control and mapping** : quality control and mapping performed by fastp and salmon.
- **Count matrix Generation**: Performed using tx2gene in R environment.
- **Data Processing**: Preprocessing and normalization of scRNA-seq data.
- **Clustering**: Identification of cell clusters using t-SNE and Seurat.
- **Pseudotime Analysis**: Inference of cell trajectories and pseudotime using Slingshot.
- **Visualization**: Generation of t-SNE plots and pseudotime trajectory plots.

## Installation

To run the analysis, you need to have R and the following packages installed:


