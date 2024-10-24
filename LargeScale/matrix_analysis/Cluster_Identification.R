library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(gridExtra)
library(grid)

islet_intr <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/integrated_data.rds")
output_fig <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/"

#---------------------------------View Clusters--------------------------------
p1 <- DimPlot(islet_intr, reduction = "umap", label = TRUE)
p2 <- DimPlot(islet_intr, reduction = 'umap', group.by = "patient_name")
p3 <- DimPlot(islet_intr, reduction = 'umap', group.by = "condition")

grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

#----------------------------Finding INS---------------------------------------
# Plotting INS expression on UMAP
FeaturePlot(islet_intr, features = "INS", min.cutoff = "q50", max.cutoff = "q90")
#Acinar gnes?

output_fig_acinar <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/Acinar_Associated _Genes/"
FeaturePlot(islet_intr, features = "AMY2B",  pt.size = 2)
FeaturePlot(islet_intr, features = "AMY2A",  pt.size = 2)
FeaturePlot(islet_intr, features = "GNMT",  pt.size = 2)
FeaturePlot(islet_intr, features = "PLIN4", pt.size = 2)
FeaturePlot(islet_intr, features = "CLPS", pt.size = 2)
FeaturePlot(islet_intr, features = "CELA3A", pt.size = 2)



output_fig_acinar <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/Acinar_Associated_Genes/"

# Plot and save AMY2B
p1 <- FeaturePlot(islet_intr, features = "AMY2B")
ggsave(filename = paste0(output_fig_acinar, "AMY2B_featureplot.png"), plot = p1)

# Plot and save AMY2A
p2 <- FeaturePlot(islet_intr, features = "AMY2A")
ggsave(filename = paste0(output_fig_acinar, "AMY2A_featureplot.png"), plot = p2)

# Plot and save GNMT
p3 <- FeaturePlot(islet_intr, features = "GNMT")
ggsave(filename = paste0(output_fig_acinar, "GNMT_featureplot.png"), plot = p3)

# Plot and save PLIN4
p4 <- FeaturePlot(islet_intr, features = "PLIN4")
ggsave(filename = paste0(output_fig_acinar, "PLIN4_featureplot.png"), plot = p4)

# Plot and save CLPS
p5 <- FeaturePlot(islet_intr, features = "CLPS")
ggsave(filename = paste0(output_fig_acinar, "CLPS_featureplot.png"), plot = p5)

# Plot and save CELA3A
p6 <- FeaturePlot(islet_intr, features = "CELA3A")
ggsave(filename = paste0(output_fig_acinar, "CELA3A_featureplot.png"), plot = p6)




#-------------------------------Identifying clusters -------------------------

DefaultAssay(islet_intr) <- "RNA"
islet_intr <- JoinLayers(islet_intr)

# # cluster 0
# markers_cluster_0 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 0,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_0) #alpha cells

islet_intr <- RenameIdents(islet_intr, '0' = 'Alpha')

# # cluster 1
# markers_cluster_1 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 1,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_1) #beta cells
islet_intr <- RenameIdents(islet_intr, '1' = 'Beta')

# cluster 2
# markers_cluster_2 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 2,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_2) #acinar cells
islet_intr <- RenameIdents(islet_intr, '2' = 'Acinar')

# cluster 3
# markers_cluster_3 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 3,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_3) #epithelial cells
islet_intr <- RenameIdents(islet_intr, '3' = 'Epithelial')

# cluster 4
# markers_cluster_4 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 4,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_4) #fibroblast or endothelial cells
islet_intr <- RenameIdents(islet_intr, '4' = 'Endothelial')
# cluster 5
# markers_cluster_5 <- FindConservedMarkers(islet_intr,
#                                           ident.1 = 5,
#                                           grouping.var = 'condition')
# 
# head(markers_cluster_5) # delta cells
islet_intr <- RenameIdents(islet_intr, '5' = 'Delta')

saveRDS(islet_intr, "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/intr_labeled.rds")

# Saving cluster labels figure
cluster_labels <- DimPlot(islet_intr, reduction = 'umap', label = TRUE)
# Save the plot to a file
ggsave(paste0(output_fig,"cluster_labels.png"), 
       plot = cluster_labels, 
       width = 10, 
       height = 8, 
       dpi = 300)

ins_plot <- FeaturePlot(islet_intr, features = "INS", min.cutoff = "q50", max.cutoff = "q90")
cluster_insulin <- plot_grid(cluster_labels, ins_plot, labels = "AUTO")

# ggsave(paste0(output_fig,"cluster_insulin.png"), 
#        plot = cluster_insulin, 
#        width = 10, 
#        height = 8, 
#        dpi = 300)

#-------------------Comparing T2D beta cells and Healthy------------------------
islet_intr$celltype <- Idents(islet_intr)

islet_intr$celltype.cnd <- paste0(islet_intr$celltype,'_', islet_intr$condition)

#relabeling based on conditions
Idents(islet_intr) <- islet_intr$celltype.cnd

# Find markers between beta cells
diseased_betacells <- FindMarkers(islet_intr, ident.1 = 'Beta_T2D', 
                                  ident.2 = 'Beta_Healthy')

head(diseased_betacells, n = 10)

saveRDS(islet_intr, "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/intr_labeled.rds")

# Finding markers between acinar cells
diseased_acinarcells <- FindMarkers(islet_intr, ident.1 = 'Acinar_T2D', 
                                  ident.2 = 'Acinar_Healthy')

head(diseased_acinarcells, n = 50)





