library(cowplot)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(fda)


BiocManager::install("monocle3")
beta_cells <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/intr_labeled.rds")
# Set up
output_fig <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/"

beta_cells <- subset(beta_cells, idents = c("Beta_T2D", "Beta_Healthy"))

Idents(beta_cells)

#---------------------------Pre processing-------------------------------------
beta_cells <- NormalizeData(beta_cells)
beta_cells <- FindVariableFeatures(beta_cells)
beta_cells <- ScaleData(beta_cells)
beta_cells <- RunPCA(beta_cells)
ElbowPlot(beta_cells)
beta_cells <- FindNeighbors(beta_cells, dims = 1:20)
beta_cells <- FindClusters(beta_cells, resolution = 0.6)
beta_cells <- RunUMAP(beta_cells, dims = 1:20, resolution = 0.6)

#---------------------------Testing------------------------------------------
# 
# beta_cells <- RunUMAP(beta_cells, dims = 1:20, resolution = c(0.03, 0.05, 0.07, 0.09))
# # beta_cells <- FindClusters(beta_cells, resolution = c(0.05, 0.06, 0.07, 0.08, 0.1))
# # 
# # DimPlot(beta_cells, group.by = "RNA_snn_res.0.1", label = TRUE)
# # DimPlot(beta_cells, group.by = "RNA_snn_res.0.08", label = TRUE)
# # DimPlot(beta_cells, group.by = "RNA_snn_res.0.07", label = TRUE)
# # DimPlot(beta_cells, group.by = "RNA_snn_res.0.06", label = TRUE)
# # DimPlot(beta_cells, group.by = "RNA_snn_res.0.05", label = TRUE)
# 
# 
# 
DimPlot(beta_cells, reduction = 'umap')
DimPlot(beta_cells, reduction = 'umap', group.by = 'condition')
# 
# 
# #------------------------------More Testing------------------------------------
# 
# # List of resolutions
# resolutions <- c(0.2, 0.4, 0.6, 0.8)
# 
# # Loop through each resolution, compute UMAP, and plot
# plots <- list()
# for (res in resolutions) {
#   # Update resolution in clustering
#   beta_cells <- FindClusters(beta_cells, resolution = res, algorithm = 3)  # algorithm 3 is Louvain
#   
#   # Run UMAP
#   beta_cells <- RunUMAP(beta_cells, dims = 1:20, 
#                         reduction.key = paste0("UMAP", res, "_"), 
#                         reduction.name = paste0("umap", res))
#   
#   # Store the plot in the list
#   plots[[paste0("Resolution ", res)]] <- DimPlot(beta_cells, 
#                                                  reduction = paste0("umap", res), 
#                                                  group.by = "seurat_clusters") +
#     ggtitle(paste("UMAP at resolution", res))
# }
# 
# # Plot all UMAPs
# 
# plot_grid(plotlist = plots, ncol = 2)


#---------------------------Monocle------------------------------------------

cds <- as.cell_data_set(beta_cells)


# # to get cell metadata
# colData(cds)
# 
# # to get gene metadata
# fData(cds)
# rownames(fData(cds))[1:10]

fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
# 
# counts(cds)



# Need monocle to use clustering information used in seurat
# Needs to be added to the cds dataset
# Monocle3 has function to clusters cells where it no only determines the clusters
# but also partitions, these are nothing but super clusters

# To do: 
# 1. add the partition information
# 2. add the clustering information from serurat
# 3. UMAP cell embeddings - coordinates


# 1. assigning the partitions

recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition


# 2. Assign cluster info

list_cluster <- beta_cells@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# 3. UMAP embeddings

cds@int_colData@listData$reducedDims$UMAP <- beta_cells@reductions$umap@cell.embeddings


#----------------------------visualizing---------------------------------------

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',  # Color cells based on cluster assignment
                                        show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
                                        cell_size = 0.8,
                                        group_label_size = 5) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )

# Print the plot
print(cluster.before.trajectory)

# This doesn't work for me, since I am not renaming the beta cells
# cluster.names <- plot_cells(cds,
#            color_cells_by = "redefined_cluster",
#            label_groups_by_cluster = FALSE,
#            group_label_size = 5)+
#   theme(legend.position = 'right')




cluster.before.trajectory

#----------------------build trajectory---------------------------------------
# use_partition = FALSE - This means that it will learn disjointed trajectories
# FALSE means that i'll do a single trajectory for everything

# cds <- learn_graph(cds, use_partition = FALSE)
cds <- learn_graph(cds)

trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_branch_points = FALSE
           
           )


trajectory_plot <- plot_cells(cds,
                              color_cells_by = 'condition',
                              label_cell_groups = FALSE,
                              cell_size = 0.8,
                              group_label_size = 50,
                              label_roots = FALSE,
                              label_leaves = TRUE) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  ) +
  labs(title = "Cell Trajectory by Condition",
       x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles


# Print the plot
print(trajectory_plot)
cluster.before.trajectory | trajectory_plot


# Saving figure
combined_plot <- plot_grid(cluster.before.trajectory, trajectory_plot, labels = "AUTO")

ggsave(paste0(output_fig, "combined_trajectory_plots.png"), 
       plot = combined_plot, 
       width = 16, 
       height = 8, 
       dpi = 300)



#----------------------------Ordering Cells------------------------------------'

# I am assuming that cluster 0 contains the root cell population, because it
# is where most of the healthy cells are clustered

cds2 <- cds
cds2 <- order_cells(cds2, reduction_method = 'UMAP', 
                   root_cells = colnames(cds2[,clusters(cds2) == 0]))


plot_cells(cds2, 
           color_cells_by = 'pseudotime',
           cell_size = 0.8,
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE)

pseudotime_plot <- plot_cells(cds2,
                              color_cells_by = 'pseudotime',
                              label_cell_groups = FALSE,
                              cell_size = 0.8,# Adjust this as needed
                              group_label_size = 50) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Cell Trajectory by Condition",
       x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles

# Print the plot
print(trajectory_plot) | pseudotime_plot

# Saving figure
pseudotime_plots <- plot_grid(cluster.before.trajectory, pseudotime_plot, labels = "AUTO")

ggsave(paste0(output_fig, "pseudotime_plots.png"), 
       plot = pseudotime_plots, 
       width = 16, 
       height = 8, 
       dpi = 300)


#---------------------------Visualization & save------------------------------

# First plot - 
#################
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',  # Color cells based on cluster assignment
                                        show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
                                        cell_size = 0.8,
                                        group_label_size = 5) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )

trajectory_plot <- plot_cells(cds,
                              color_cells_by = 'condition',
                              label_cell_groups = FALSE,
                              cell_size = 0.8,
                              group_label_size = 50,
                              label_roots = FALSE,
                              label_leaves = TRUE) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  ) +
  labs(title = "Cell Trajectory by Condition",
       x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles

# Saving figure
combined_plot <- plot_grid(cluster.before.trajectory, trajectory_plot, labels = "AUTO")

ggsave(paste0(output_fig, "combined_trajectory_plots.png"), 
       plot = combined_plot, 
       width = 16, 
       height = 8, 
       dpi = 300)

# Pseudotime plots healhty vs t2d and pseudotime
########################

# Left - healthy vs t2d

plot.conditions<- plot_cells(cds,
                             color_cells_by = 'condition',  # Color cells based on cluster assignment
                             label_cell_groups = FALSE,
                             show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
                             cell_size = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )


# Pseudotime plot
pseudotime_plot <- plot_cells(cds2,
                              color_cells_by = 'pseudotime',
                              label_cell_groups = FALSE,
                              cell_size = 0.8,
                              label_groups_by_cluster = FALSE,
                              # label_branch_points = FALSE,
                              # label_roots = FALSE
                              ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )+
  labs(title = "Cell Trajectory by Condition",
       x = "UMAP 1", y = "UMAP 2")

# Saving figure
pseudotime_plots <- plot_grid(plot.conditions, pseudotime_plot, labels = "AUTO")

ggsave(paste0(output_fig, "pseudotime_plots.png"), 
       plot = pseudotime_plots, 
       width = 16, 
       height = 8, 
       dpi = 300)


#----------------------Capturing age data -------------------------------------

library(readxl)

metadata <- read_excel("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/metadata/PancDB_Donors.xlsx")
metadata$SampleAge <- as.numeric(gsub("yo", "", metadata$SampleAge))



cds_col_data <- colData(cds)
dim(cds_col_data)

unique(cds_col_data$patient_name)
unique(metadata$DonorID)
summary(metadata$SampleAge)
head(metadata)

# Create a new 'age' column in cds@colData and initialize it with NA
cds@colData$age <- NA

# Iterate through each row in cds@colData
for (i in seq_len(nrow(cds@colData))) {
  # Find the matching donor ID in the metadata
  donor_match <- metadata$DonorID == cds@colData$patient_name[i]
  
  # If there's a match, update the 'age' field for that row in cds@colData
  if (any(donor_match)) {
    cds@colData$age[i] <- metadata$SampleAge[donor_match]
  }
}


Age <- plot_cells(cds,
           color_cells_by = 'age',  # Color cells based on cluster assignment
           label_cell_groups = FALSE,
           # show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
           cell_size = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )


# Saving figure
AgenPseudotime_plots <- plot_grid(Age, pseudotime_plot, labels = "AUTO")

ggsave(paste0(output_fig, "AgenPseudotime_plots.png"), 
       plot = AgenPseudotime_plots, 
       width = 16, 
       height = 8, 
       dpi = 300)


#-------------------------Age related analysis---------------------------------

sort(unique(cds@colData$age))




