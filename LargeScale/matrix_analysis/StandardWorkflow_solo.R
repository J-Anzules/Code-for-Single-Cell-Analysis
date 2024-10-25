library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(SeuratData)

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

FeatureScatter(T2D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method='lm') +
  geom_hline(yintercept = 1850)+
  geom_vline(xintercept = 300000)

FeatureScatter(healthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method='lm')


#--------------------------Standard Processing -------------------------------


# Normalizing
T2D <- NormalizeData(T2D)
healthy <- NormalizeData(healthy)

# Highly variable genes
T2D <- FindVariableFeatures(T2D, selection.method = "vst", nfeatures = 2000)
healthy <- FindVariableFeatures(healthy, selection.method = "vst", nfeatures = 2000)

# Top 10 of both
top10_T2D <- head(VariableFeatures(T2D), 10)
top10_healthy <- head(VariableFeatures(healthy), 10)

### Figures - Highly variable genes
highly_variable_t2d <- paste0(output_fig, "highly_variable_t2d.png")
highly_variable_healthy <- paste0(output_fig, "highly_variable_healthy.png")

# T2D
png(filename = highly_variable_t2d, width = 1000, height = 800)
variable_t2d <- VariableFeaturePlot(T2D)
LabelPoints(plot = variable_t2d, points = top10_T2D)
dev.off()

# healthy
png(filename = highly_variable_healthy, width = 1000, height = 800)
variable_healthy <- VariableFeaturePlot(healthy)
LabelPoints(plot = variable_healthy, points = top10_healthy)
dev.off()

#---------------------Accounting for variations-------------------------------
# Cell cycles and batch effects

all.genes_T2D <- rownames(T2D)
T2D <- ScaleData(T2D, features = all.genes_T2D)

all.genes_healthy <- rownames(healthy)
healthy <- ScaleData(healthy, features = all.genes_healthy)

#--------------------Linear Dimensionality reduction PCA------------------------
#Identify the sources of heterogeneity in the dataset

T2D <- RunPCA(T2D, features = VariableFeatures(object = T2D))
healthy <- RunPCA(healthy, features = VariableFeatures(object = healthy))

# dims = first 5 PCA but only the first 5 (nfeatures) features
print(T2D[["pca"]], dims=1:5, nfeatures = 5)
DimHeatmap(T2D, dims = 1, cells = 500, balanced = TRUE)

# Elbow plot
# Choosing only those statistically significant components that captures the
# majority of the signal in our downstream analysis
#Ranked by the percentage of the variance explained

ElbowPlot(T2D) #12
ElbowPlot(healthy) #12

# These could adversely effect the downstream analysis


#----------------------Clustering based on PCA----------------------------------
T2D <- FindNeighbors(T2D, dims = 1:12)
healthy <- FindNeighbors(healthy, dims = 1:12)

# Checking out the resolutions
T2D <- FindClusters(T2D, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(T2D@meta.data$RNA_snn_res.0.1)
healthy <- FindClusters(healthy, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

DimPlot(T2D, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(healthy, group.by = "RNA_snn_res.0.7", label = TRUE)

# Setting identify of clusters
Idents(T2D) <- "RNA_snn_res.0.5"
Idents(healthy) <- "RNA_snn_res.0.7"

#--------------------------------UMAP-------------------------------------------
# reticulate::py_install(packages = 'umap-learn')
library(umap)

T2D <- RunUMAP(T2D, dims = 1:12)
healthy <- RunUMAP(healthy, dims = 1:12)

DimPlot(T2D, reduction = "umap")
DimPlot(T2D, reduction = "umap", group.by = "patient_name")

DimPlot(healthy, reduction = "umap")
DimPlot(healthy, reduction = "umap", group.by = "patient_name")

unique(T2D@meta.data$patient_name)


#-----------------------------Saving-------------------------------------------
SaveH5Seurat(object = T2D,
              filename = "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/T2D_clean_clustered") 


#--------------------------Testing Merging -------------------------------------

# Merge Seurat objects
merged_object <- merge(seurat_object1, y = seurat_object2, 
                       add.cell.ids = c("Set1", "Set2"), 
                       project = "Merged_Dataset")
