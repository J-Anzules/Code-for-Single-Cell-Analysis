#---
#  title: "Profiling Beta cell heterogeneity in T2D"
#author: "Jonathan Anzules"
#date: "2023-12-29"
#output: html_document
#---
  

library(Seurat)
library(dplyr)
library(ggplot2)
library(slingshot) # Pseudotime analysis
library(grDevices) # Pseudotime colors
library(RColorBrewer) # Pseudotime colors
library(monocle3) # Pseudotime 2
library(SeuratWrappers) # For cell_data_set

HQ_counts <- read.csv("/Documents and Settings/jonan/Documents/1Work/scWork/Beta Cell Study/Data/FirstCountMatrix.csv", row.names = 1)

# Creating a Seurat object
pancr_cells <- CreateSeuratObject(counts = HQ_counts)

# Normalizing to 10000 count
pancr_cells <- NormalizeData(pancr_cells, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling
all.genes <- rownames(pancr_cells)
pancr_cells <- ScaleData(pancr_cells, features = all.genes)

# Identifying variable features
pancr_cells <- FindVariableFeatures(pancr_cells, selection.method = "vst", nfeatures = 2000)

pancr_cells <- RunPCA(pancr_cells, features = VariableFeatures(object = pancr_cells))

pancr_cells <- FindNeighbors(pancr_cells, dims = 1:10)
pancr_cells <- FindClusters(pancr_cells, resolution = 0.4)

beta_cells <- subset(pancr_cells, idents = c("2", "4", "5"))

# Save the new dataframe for beta cells
saveRDS(beta_cells, file = "/Documents and Settings/jonan/Documents/1Work/scWork/Beta Cell Study/Data/beta_cells.rds")




# Convert Seurat object to Monocle's CellDataSet format
cds <- as.cell_data_set(beta_cells)

# Preprocess the data (normalization and feature selection)
cds <- preprocess_cds(cds, num_dim = 50)

# Reduce the dimensionality using UMAP (or PCA, if preferred)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Cluster the cells in Monocle
cds <- cluster_cells(cds)

# Learn the trajectory graph
cds <- learn_graph(cds)

# Order cells in pseudotime
cds <- order_cells(cds)

# Plot the pseudotime trajectory
plot_cells(cds, color_cells_by = "pseudotime")



