library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(gridExtra)
library(grid)

T2D <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_T2D.rds")
healthy <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/count_matrix_healthy.rds")
T2D$condition <- 'T2D'
healthy$condition <- 'Healthy'

#------------------------Quality Control n Merging------------------------------
T2D[["percent.mt"]] <- PercentageFeatureSet(T2D, pattern = "^MT-")
healthy[["percent.mt"]] <- PercentageFeatureSet(healthy, pattern = "^MT-")

T2D <- subset(T2D, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000 & nCount_RNA < 1850000 & percent.mt < 78)
# dim(T2D) Original - 42148  3862
# dim(T2D) New      - 42148  3435
healthy <- subset(healthy, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & nCount_RNA < 1850000 & percent.mt < 78)
# dim(healthy) Original - 42148  7341
# dim(healthy) New      - 42148  5868

# Merge the two Seurat objects with an additional layer of identifiers
Islets_merged <- merge(T2D, y = healthy, add.cell.ids = c("T2D", "Healthy"), 
                       project = "Combined_Islets")

#------------------------Standard Workflow-------------------------------------

Islets_merged <- NormalizeData(Islets_merged)
Islets_merged <- FindVariableFeatures(Islets_merged)
Islets_merged <- ScaleData(Islets_merged)
Islets_merged <- RunPCA(Islets_merged)
ElbowPlot(Islets_merged)
Islets_merged <- FindNeighbors(Islets_merged, dims = 1:20)
Islets_merged <- FindClusters(Islets_merged)
Islets_merged <- RunUMAP(Islets_merged, dims = 1:20)



#-------------------------Visualizing batch effect------------------------------

head(Islets_merged@meta.data)

p1 <- DimPlot(Islets_merged, reduction = "umap")
p2 <- DimPlot(Islets_merged, reduction = 'umap', group.by = "patient_name")
p3 <- DimPlot(Islets_merged, reduction = 'umap', group.by = "condition")

combined_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

png("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/before_integration.png",
    width = 2000, height = 1200)

grid.draw(combined_plot)
dev.off()

#Signs of technical variation causes clustering, not biological variation



#--------------------Performing integration------------------------------------
obj.list <- SplitObject(Islets_merged, split.by = 'patient_name')

for( i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# Selecting Integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# Find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                       anchor.features = features, k.filter = )

# Integrate data
islet_integrated <- IntegrateData(anchorset = anchors)


#-------------------Standard Workflow -----------------------------------------
islet_integrated <- ScaleData(object = islet_integrated)
islet_integrated <- RunPCA(object = islet_integrated)
ElbowPlot(islet_integrated)
islet_integrated <- FindNeighbors(islet_integrated, dims = 1:11)

# Checking out the resolutions
islet_integrated <- FindClusters(islet_integrated, resolution = c(0.05, 0.06, 0.07, 0.08))
View(islet_integrated@meta.data)

DimPlot(islet_integrated, group.by = "integrated_snn_res.0.05", label = TRUE)

islet_integrated <- RunUMAP(islet_integrated, dims = 1:7, resolution = 0.05)

# Saving Integrated data
saveRDS(islet_integrated, file = "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/integrated_data.rds")


#-------------------------Visualizing-----------------------------------------


p1 <- DimPlot(islet_integrated, reduction = "umap", label = TRUE)
p2 <- DimPlot(islet_integrated, reduction = 'umap', group.by = "patient_name")
p3 <- DimPlot(islet_integrated, reduction = 'umap', group.by = "condition")

grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

combined_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

png("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/after_integration.png",
    width = 2000, height = 1200)

grid.draw(combined_plot)
dev.off()

DefaultAssay(islet_integrated)











