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


# Generate Metadata
#####

# Define the base directory
base_dir <- "C://Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/quant_files/upenn1/"

# List all quant.sf files
quant_files <- list.files(base_dir, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)


###################################################################################
###################################################################################
###################################################################################
###################################################################################

# I need to add an if statment to extract the meta data for this file type now:
# HPAP-058_scRNA_SSq2-090P01_quant
# These are the first bath of scRNA-seq data that I downloaded. 

###################################################################################
###################################################################################
###################################################################################
###################################################################################

# Extract patient names and cell identifiers
# extract_info <- function(filepath) {
#   # Extract patient name and cell ID using regex
#   match <- regmatches(filepath, regexec("([^/]+)_scRNA_(\\d+)_quant/quant.sf", filepath))
#   patient_name <- match[[1]][2]
#   cell_id <- match[[1]][3]
#   return(c(filepath, patient_name, cell_id))
# }
# metadata_list <- lapply(quant_files, extract_info)

#str(metadata_list)
# Converting list of lists to dataframe
metadata_df <- do.call(rbind, metadata_list) %>% as.data.frame()
colnames(metadata_df) <- c("quant_file", "patient_name", "cell_id")

# Grabbing files address
files <- metadata_df$quant_file

#DELETE: 
files_short <- files[1:300]
metadata_df_short <- head(metadata_df, n = 300)

#Reading in the tx2gene dataframe
tx2gene <- read.csv("/Documents and Settings/jonan/Documents/Genomes/Homo-Sapiens/GRCh38/tx2gene-Mappings")

# Import the quant files using tximport
txi_short <- tximport(files_short, type = "salmon", tx2gene = tx2gene)

# Assign cell IDs from metadata as column names of the count matrix
colnames(txi_short$counts) <- metadata_df_short$cell_id

# Assign patient names as a new column in the count matrix metadata
txi_short$meta <- metadata_df_short

str(txi_short)




