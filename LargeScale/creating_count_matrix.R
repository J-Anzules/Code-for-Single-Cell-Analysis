library(tximport) # Importing quant_sf files
library(GenomicFeatures) #Useful genomics tools
library(BiocManager)
library(rhdf5) #To make an hdf5 files
library(Rtsne)
library(DESeq2) # For normalizing
library(scran) # normalization based on cell clusters
library(SingleCellExperiment) #For holding single-cell RNA_seq data
library(biomaRt) # Switching names from ensemble ID to gene symbol
library(Seurat)
library(ggplot2)
library(cowplot) # For violin plots
library(dplyr)

"C://Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/quant_files/upenn1/HPAP-001_scRNA_43690_quant/quant.sf"


# Organize quant.sf files
#####

# Define the base directory
base_dir <- "C://Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/quant_files/upenn1/"

# List all quant.sf files
quant_files <- list.files(base_dir, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# Extract patient names and cell identifiers
extract_info <- function(filepath) {
  # Extract patient name and cell ID using regex
  match <- regmatches(filepath, regexec("([^/]+)_scRNA_(\\d+)_quant/quant.sf", filepath))
  patient_name <- match[[1]][2]
  cell_id <- match[[1]][3]
  return(c(filepath, patient_name, cell_id))
}

quant_files[1]
a <- regmatches(quant_files[1], regexec("([^/]+)_scRNA_(\\d+)_quant/quant.sf", quant_files[1]))
a[[1]]




