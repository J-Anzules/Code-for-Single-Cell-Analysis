library(Seurat)
library(ggplot2)

# library(GenomicFeatures) #Useful genomics tools
# library(BiocManager)
# library(rhdf5) #To make an hdf5 files
# library(Rtsne)
# library(DESeq2) # For normalizing
# library(scran) # normalization based on cell clusters
# library(SingleCellExperiment) #For holding single-cell RNA_seq data
# library(biomaRt) # Switching names from ensemble ID to gene symbol
# library(cowplot) # For violin plots


T2D <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_T2D.rds")
healthy <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_healthy.rds")

# Set up
output_fig <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/"

#----------------------------------Cleaning------------------------------------
T2D[["percent.mt"]] <- PercentageFeatureSet(T2D, pattern = "^MT-")
healthy[["percent.mt"]] <- PercentageFeatureSet(healthy, pattern = "^MT-")

qc_output_t2d <- paste0(output_fig, "QC_violin_t2d.png")
qc_output_healthy <- paste0(output_fig, "QC_violin_healthy.png")

# Visualize QC metrics as a violin plot - t2d
png(filename = qc_output_t2d, width = 1000, height = 800)
VlnPlot(T2D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Visualize QC metrics as a violin plot - healthy
png(filename = qc_output_healthy, width = 1000, height = 800)
VlnPlot(healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

T2D <- subset(T2D, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000 & nCount_RNA < 1850000 & percent.mt < 78)
# dim(T2D) Original - 42148  3862
# dim(T2D) New      - 42148  3435
healthy <- subset(healthy, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & nCount_RNA < 1850000 & percent.mt < 78)
# dim(healthy) Original - 42148  7341
# dim(healthy) New      - 42148  5868

#----------------------------------Cleaning------------------------------------









































