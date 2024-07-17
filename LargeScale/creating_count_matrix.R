library(dplyr)
library(tximport) # Importing quant_sf files
library(biomaRt) # Converting ensembleID to geneID
library(Seurat)

# library(GenomicFeatures) #Useful genomics tools
# library(BiocManager)
# library(rhdf5) #To make an hdf5 files
# library(Rtsne)
# library(DESeq2) # For normalizing
# library(scran) # normalization based on cell clusters
# library(SingleCellExperiment) #For holding single-cell RNA_seq data
# library(biomaRt) # Switching names from ensemble ID to gene symbol

# library(ggplot2)
# library(cowplot) # For violin plots



# Generate Metadata
#####

# Define the base directory
# base_dir <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/quant_files_T2D/upenn1/"
base_dir <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/quant_files_healthy/upenn1/"

# List all quant.sf files
quant_files <- list.files(base_dir, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# quant_files <- head(quant_files, n = 10)
extract_info <- function(filepath) {
  # Extract patient name and cell ID using regex for the first pattern
  match1 <- regmatches(filepath, regexec("([^/]+)_scRNA_(\\d+)_quant/quant.sf", filepath))
  
  # Extract patient name and cell ID using regex for the second pattern
  match2 <- regmatches(filepath, regexec("([^/]+)_scRNA_SSq2-([^_]+)_quant", filepath))
  
  if (length(match1[[1]]) > 1) {
    patient_name <- match1[[1]][2]
    cell_id <- match1[[1]][3]
  } else if (length(match2[[1]]) > 1) {
    patient_name <- match2[[1]][2]
    cell_id <- match2[[1]][3]
  } else {
    patient_name <- NA
    cell_id <- NA
  }
  
  return(c(filepath, patient_name, cell_id))
}

metadata_list <- lapply(quant_files, extract_info)


#str(metadata_list)
# Converting list of lists to dataframe
metadata_df <- do.call(rbind, metadata_list) %>% as.data.frame()
colnames(metadata_df) <- c("quant_file", "patient_name", "cell_id")

# Grabbing files address
files <- metadata_df$quant_file

#Reading in the tx2gene dataframe
tx2gene <- read.csv("/Documents and Settings/jonan/Documents/Genomes/Homo-Sapiens/GRCh38/tx2gene-Mappings")

# Import the quant files using tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


#--------------- Replacing Ensemble ID's with gene names ----------------------#
########################
ensemble_ids <- gsub("\\.\\d+$", "", rownames(txi$counts))
rownames(txi$counts) <- gsub("\\.\\d+$", "", rownames(txi$counts))

# Connecting to Ensembl and retrieve gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

converted_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = ensemble_ids,
                       mart = ensembl)


# Create a mapping from Ensembl IDs to Gene Names
id_to_name <- setNames(converted_ids$external_gene_name, converted_ids$ensembl_gene_id)

# Replace row names in norm_counts with gene names, using Ensembl IDs as keys
txi_rownames <- rownames(txi$counts)
rownames(txi$counts) <- ifelse(txi_rownames %in% names(id_to_name), id_to_name[txi_rownames], NA)

# Filter out rows where row names are NA or empty
txi$counts <- txi$counts[!is.na(rownames(txi$counts)) & rownames(txi$counts) != "", ] #878 NA's

counts <- txi$counts
colnames(counts) <- metadata_df$cell_id
rownames(counts) <- make.unique(rownames(counts), sep="-")

#----------------------------------------------------------------------------
#                   Making Seurat Object

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

# Assign patient names as a new column in the count matrix metadata
metadata_df$quant_file <- NULL
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_df)

saveRDS(seurat_obj, file = "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_T2D.rds")
saveRDS(seurat_obj, file = "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_healthy.rds")