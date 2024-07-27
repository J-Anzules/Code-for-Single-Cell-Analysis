library(cowplot)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(fda)


integrated_data <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/intr_labeled.rds")
# Set up
output_fig_ins <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/ins_cells/"

#----------------------------Subsetting insulin cells --------------------------

qtle <- "q50"
FeaturePlot(integrated_data, features = "INS", min.cutoff = qtle, max.cutoff = "q90")+
  ggtitle(paste("Expression of INS", qtle)) + 
  theme(plot.title = element_text(hjust = 0.5))


ins_expression <- FetchData(integrated_data, vars = "INS")

# Calculate the 40th and 90th quantiles of INS expression
quantile_40 <- quantile(ins_expression$INS, 0.4)
quantile_90 <- quantile(ins_expression$INS, 0.9)

# Create a logical vector that indicates which cells fall within the quantile range
cells_in_range <- (ins_expression$INS > quantile_40 & ins_expression$INS < quantile_90)

# length(cells_in_range)


# Get the cell names (column names) that correspond to TRUE in the logical vector
cell_names_to_keep <- colnames(integrated_data)[cells_in_range]

# Subset the Seurat object by these cell names
ins_cells <- integrated_data[, cell_names_to_keep]
# saveRDS(ins_cells, "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/insulin_cells.rds")
#---------------------------Pre processing-------------------------------------
ins_cells <- NormalizeData(ins_cells)
ins_cells <- FindVariableFeatures(ins_cells)
ins_cells <- ScaleData(ins_cells)
ins_cells <- RunPCA(ins_cells)
# ElbowPlot(ins_cells)
ins_cells <- FindNeighbors(ins_cells, dims = 1:10)
ins_cells <- FindClusters(ins_cells, resolution = 0.6)
ins_cells <- RunUMAP(ins_cells, dims = 1:10, resolution = 0.6)

#---------------------------Testing------------------------------------------
# 
# ins_cells <- RunUMAP(ins_cells, dims = 1:20, resolution = c(0.03, 0.05, 0.07, 0.09))
# ins_cells <- FindClusters(ins_cells, resolution = c(0.05, 0.06, 0.07, 0.08, 0.1))
# 
# DimPlot(ins_cells, group.by = "RNA_snn_res.0.1", label = TRUE)
# DimPlot(ins_cells, group.by = "RNA_snn_res.0.08", label = TRUE)
# DimPlot(ins_cells, group.by = "RNA_snn_res.0.07", label = TRUE)
# DimPlot(ins_cells, group.by = "RNA_snn_res.0.06", label = TRUE)
# DimPlot(ins_cells, group.by = "RNA_snn_res.0.05", label = TRUE)
# 
# 
# 
# DimPlot(ins_cells, reduction = 'umap')
# DimPlot(ins_cells, reduction = 'umap', group.by = 'condition')


# #------------------------------More Testing------------------------------------
# 
# # List of resolutions
# resolutions <- c(0.2, 0.4, 0.6, 0.8)
# 
# # Loop through each resolution, compute UMAP, and plot
# plots <- list()
# for (res in resolutions) {
#   # Update resolution in clustering
#   ins_cells <- FindClusters(ins_cells, resolution = res, algorithm = 3)  # algorithm 3 is Louvain
# 
#   # Run UMAP
#   ins_cells <- RunUMAP(ins_cells, dims = 1:20,
#                         reduction.key = paste0("UMAP", res, "_"),
#                         reduction.name = paste0("umap", res))
# 
#   # Store the plot in the list
#   plots[[paste0("Resolution ", res)]] <- DimPlot(ins_cells,
#                                                  reduction = paste0("umap", res),
#                                                  group.by = "seurat_clusters") +
#     ggtitle(paste("UMAP at resolution", res))
# }
# 
# # Plot all UMAPs
# 
# plot_grid(plotlist = plots, ncol = 2)


#---------------------------Monocle------------------------------------------

cds_ins <- as.cell_data_set(ins_cells)


# # to get cell metadata
# colData(cds_ins)
# 
# # to get gene metadata
# fData(cds_ins)
# rownames(fData(cds_ins))[1:10]

fData(cds_ins)$gene_short_name <- rownames(fData(cds_ins))

# to get counts
# 
# counts(cds_ins)



# Need monocle to use clustering information used in seurat
# Needs to be added to the cds_ins dataset
# Monocle3 has function to clusters cells where it no only determines the clusters
# but also partitions, these are nothing but super clusters

# To do: 
# 1. add the partition information
# 2. add the clustering information from serurat
# 3. UMAP cell embeddings - coordinates


# 1. assigning the partitions

recreate.partition <- c(rep(1, length(cds_ins@colData@rownames)))
names(recreate.partition) <- cds_ins@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_ins@clusters$UMAP$partitions <- recreate.partition


# 2. Assign cluster info

list_cluster <- ins_cells@active.ident
cds_ins@clusters$UMAP$clusters <- list_cluster

# 3. UMAP embeddings

cds_ins@int_colData@listData$reducedDims$UMAP <- ins_cells@reductions$umap@cell.embeddings


#----------------------------visualizing---------------------------------------

cluster.before.trajectory <- plot_cells(cds_ins,
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
# cluster.names <- plot_cells(cds_ins,
#            color_cells_by = "redefined_cluster",
#            label_groups_by_cluster = FALSE,
#            group_label_size = 5)+
#   theme(legend.position = 'right')




cluster.before.trajectory

#----------------------build trajectory---------------------------------------
# use_partition = FALSE - This means that it will learn disjointed trajectories
# FALSE means that i'll do a single trajectory for everything

# cds_ins <- learn_graph(cds_ins, use_partition = FALSE)
cds_ins <- learn_graph(cds_ins)

# trajectory <- plot_cells(cds_ins,
#            color_cells_by = 'cluster',
#            label_branch_points = FALSE
#            
#            )
# 
# 
# trajectory_plot <- plot_cells(cds_ins,
#                               color_cells_by = 'condition',
#                               label_cell_groups = FALSE,
#                               cell_size = 0.8,
#                               group_label_size = 50,
#                               label_roots = FALSE,
#                               label_leaves = TRUE) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   ) +
#   labs(title = "Cell Trajectory by Condition",
#        x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles
# 
# 
# # Print the plot
# print(trajectory_plot)
# cluster.before.trajectory | trajectory_plot
# 
# 
# # Saving figure
# combined_plot <- plot_grid(cluster.before.trajectory, trajectory_plot, labels = "AUTO")
# 
# ggsave(paste0(output_fig_ins, "combined_trajectory_plots.png"), 
#        plot = combined_plot, 
#        width = 16, 
#        height = 8, 
#        dpi = 300)



#----------------------------Ordering Cells------------------------------------'

# I am assuming that cluster 0 contains the root cell population, because it
# is where most of the healthy cells are clustered

cds_ins2 <- cds_ins
cds_ins2 <- order_cells(cds_ins2, reduction_method = 'UMAP', 
                   root_cells = colnames(cds_ins2[,clusters(cds_ins2) == 5]))


# plot_cells(cds_ins2, 
#            color_cells_by = 'pseudotime',
#            cell_size = 0.8,
#            label_groups_by_cluster = FALSE,
#            label_branch_points = FALSE,
#            label_roots = FALSE)
# 
# pseudotime_plot <- plot_cells(cds_ins2,
#                               color_cells_by = 'pseudotime',
#                               label_cell_groups = FALSE,
#                               cell_size = 0.8,# Adjust this as needed
#                               group_label_size = 50) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   )+
#   labs(title = "Cell Trajectory by Condition",
#        x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles
# 
# # Print the plot
# # print(trajectory_plot) | pseudotime_plot
# 
# # Saving figure
# pseudotime_plots <- plot_grid(cluster.before.trajectory, pseudotime_plot, labels = "AUTO")
# 
# ggsave(paste0(output_fig_ins, "pseudotime_plots.png"), 
#        plot = pseudotime_plots, 
#        width = 16, 
#        height = 8, 
#        dpi = 300)


#---------------------------Visualization & save------------------------------

# First plot - 
#################
# cluster.before.trajectory <- plot_cells(cds_ins,
#                                         color_cells_by = 'cluster',  # Color cells based on cluster assignment
#                                         show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
#                                         cell_size = 0.8,
#                                         group_label_size = 5) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   )
# 
# trajectory_plot <- plot_cells(cds_ins,
#                               color_cells_by = 'condition',
#                               label_cell_groups = FALSE,
#                               cell_size = 0.8,
#                               group_label_size = 50,
#                               label_roots = FALSE,
#                               label_leaves = TRUE) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   ) +
#   labs(title = "Cell Trajectory by Condition",
#        x = "UMAP 1", y = "UMAP 2")  # Customize your labels and titles
# 
# # Saving figure
# combined_plot <- plot_grid(cluster.before.trajectory, trajectory_plot, labels = "AUTO")
# 
# ggsave(paste0(output_fig_ins, "combined_trajectory_plots.png"), 
#        plot = combined_plot, 
#        width = 16, 
#        height = 8, 
#        dpi = 300)

# Pseudotime plots healhty vs t2d and pseudotime
########################

# Left - healthy vs t2d

# plot.conditions<- plot_cells(cds_ins,
#                              color_cells_by = 'condition',  # Color cells based on cluster assignment
#                              label_cell_groups = FALSE,
#                              show_trajectory_graph = FALSE,  # Ensure trajectories are not shown
#                              cell_size = 0.8) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   )
# 
# 
# # Pseudotime plot
# pseudotime_plot <- plot_cells(cds_ins2,
#                               color_cells_by = 'pseudotime',
#                               label_cell_groups = FALSE,
#                               cell_size = 0.8,
#                               label_groups_by_cluster = FALSE,
#                               # label_branch_points = FALSE,
#                               # label_roots = FALSE
#                               ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
#     legend.background = element_rect(fill = "white", color = NA), # White background for legend
#     legend.position = "right"
#   )+
#   labs(title = "Cell Trajectory by Condition",
#        x = "UMAP 1", y = "UMAP 2")
# 
# # Saving figure
# pseudotime_plots <- plot_grid(plot.conditions, pseudotime_plot, labels = "AUTO")
# 
# ggsave(paste0(output_fig_ins, "pseudotime_plots.png"), 
#        plot = pseudotime_plots, 
#        width = 16, 
#        height = 8, 
#        dpi = 300)


#----------------------Capturing age data -------------------------------------

library(readxl)

metadata <- read_excel("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/metadata/PancDB_Donors.xlsx")
metadata$SampleAge <- as.numeric(gsub("yo", "", metadata$SampleAge))



cds_ins_col_data <- colData(cds_ins)
dim(cds_ins_col_data)

unique(cds_ins_col_data$patient_name)
unique(metadata$DonorID)
summary(metadata$SampleAge)
head(metadata)

# Create a new 'age' column in cds_ins@colData and initialize it with NA
cds_ins@colData$age <- NA

# Iterate through each row in cds_ins@colData
for (i in seq_len(nrow(cds_ins@colData))) {
  # Find the matching donor ID in the metadata
  donor_match <- metadata$DonorID == cds_ins@colData$patient_name[i]
  
  # If there's a match, update the 'age' field for that row in cds_ins@colData
  if (any(donor_match)) {
    cds_ins@colData$age[i] <- metadata$SampleAge[donor_match]
  }
}


Age <- plot_cells(cds_ins,
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

ggsave(paste0(output_fig_ins, "AgenPseudotime_plots.png"), 
       plot = AgenPseudotime_plots, 
       width = 16, 
       height = 8, 
       dpi = 300)


Cond_Age_pseudo <- plot_grid(cluster.before.trajectory, plot.conditions, Age, pseudotime_plot, labels = "AUTO" )
ggsave(paste0(output_fig_ins, "CondAgenPseudotime_plots.png"), 
       plot = Cond_Age_pseudo, 
       width = 16, 
       height = 8, 
       dpi = 300)

# 3 year old data
###################

library(monocle3)

# Add a new column to colData to indicate if the cell is from a person aged 3
cds_ins@colData$age_group <- ifelse(cds_ins@colData$age == 3, "Age 3", "Other Ages")

# Check the new column
table(cds_ins@colData$age_group)

sort(unique(cds_ins@colData$age))
# Plot cells with color differentiation for age 3
plot <- plot_cells(cds_ins,
                   color_cells_by = 'age_group',  # Color cells based on new age group
                   label_cell_groups = FALSE,
                   cell_size = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )

# Print the plot
print(plot)

# Bins
###################
# Define the age bins
age_bins <- cut(cds_ins@colData$age,
                breaks = c(0, 5, 23, 34, 49, 70),
                right = TRUE, # right = TRUE means the intervals are closed on the right
                labels = c("0-5", "6-23", "24-34", "42-49", "50-70"),
                include.lowest = TRUE)

# Add the age bins as a new column to colData
cds_ins@colData$age_bin <- age_bins

# Check the distribution of age bins
table(cds_ins@colData$age_bin)

# Plot cells with color differentiation based on age bins
plot <- plot_cells(cds_ins,
                   color_cells_by = 'age_bin',  # Use the new age bin factor
                   label_cell_groups = FALSE,
                   cell_size = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# Print the plot
print(plot)

#-----------------------------Further ordering cells---------------------------

# Assuming 'cds' is already preprocessed and ready for testing
results <- list()
k_values <- c(5, 10, 25)  # Example range of k values

for (k in k_values) {
  cat("Testing with k =", k, "\n")
  res <- graph_test(cds_ins2, neighbor_graph = "principal_graph", 
                    reduction_method = "UMAP", k = k, verbose = TRUE,
                    cores = 10)
  results[[as.character(k)]] <- res
  # Optionally save or plot results here
}


traj_results <- graph_test(cds_ins2, neighbor_graph = "principal_graph", 
                           reduction_method = "UMAP", k = 10, verbose = TRUE,
                           cores = 10)


ordered_results <-traj_results %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head(n=32)

most <- rownames(ordered_results)[1:4]
most2 <- rownames(ordered_results)[5:8]
most3 <- rownames(ordered_results)[9:12]
#---------------------Visualizing the top 32 genes------------------------------

FeaturePlot(ins_cells, features = most)
FeaturePlot(ins_cells, features = most2)
FeaturePlot(ins_cells, features = "INS", pt.size=1.5)

genes <- rownames(ordered_results)

for (i in seq(1, length(genes), by = 4)) {
  # Get the current set of four genes
  current_genes <- genes[i:min(i + 3, length(genes))]
  
  # Generate the FeaturePlot
  plot <- FeaturePlot(ins_cells, features = current_genes)
  
  # Define the file name
  file_name <- paste0(output_fig_ins, "Top32Genes/", "FeaturePlot_", i, "_to_", min(i + 3, length(genes)), ".png")
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 10, height = 6, dpi = 300)
}


#---------------------Gene set Enrichment Analysis-----------------------------
library(dplyr)

colnames(traj_results)
dfGSEA <- traj_results %>%
  filter(status == 'OK') %>%  # filter only successful results
  arrange(desc(morans_I)) %>%  # Sort genes based on Moran's I statistic in descending order
  mutate(rank = row_number()) %>%  # Assign ranks based on Moran's I value
  dplyr::select(gene_short_name, morans_I, p_value, rank)

# Create a named vector of ranks for GSEA
ranks <- setNames(dfGSEA$rank, dfGSEA$gene_short_name)

# Check for duplicates and resolve them by retaining the highest value
duplicates <- names(ranks)[duplicated(names(ranks)) | duplicated(names(ranks), fromLast = TRUE)]
if (length(duplicates) > 0) {
  ranks <- tapply(ranks, INDEX = names(ranks), FUN = max)
}


# Import .gmt file with pathway information by gene symbol (Gene sets derived from the GO Biological Process ontology)
## Obtained from GSEA 3-6-2023 http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
go <- gmtPathways("C:/Users/jonan/Documents/Tyseq/Data/c5.go.v2023.1.Hs.symbols.gmt")

library(fgsea) 

# Run GSEA
set.seed(06212022)  # Setting a seed for reproducibility
fgseaRes <- fgsea(pathways = go,
                  stats = ranks,
                  minSize = 6,
                  maxSize = 500,
                  nproc = 1)


# Sort the results by adjusted p-value in ascending order
sorted_fgseaRes <- fgseaRes %>%
  arrange(padj) %>%  # Sorting by padj, where lower values indicate higher statistical significance
  # dplyr::select(pathway, padj) %>%
  head(n=20)


sorted_fgseaRes
head(ranks, n=50)

#---------------------------------Visualizing results--------------------------
################# Preparing Data
library(Matrix)
library(pheatmap)

# Select the top 20 pathways
topPathways <- fgseaRes %>% 
  arrange(padj) %>%
  mutate(NES_sign = ifelse(NES > 0, "Positive", "Negative")) %>%
  head(n = 20)


# Expression matrix
expression_matrix <- assay(cds_ins2, "counts") # Adjust to your specific assay type if necessary


################# BAR plot
# Create the bar plot
pathways <- ggplot(topPathways, aes(reorder(pathway, NES), NES, fill = NES_sign)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 20 Enriched Pathways") +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "firebrick"))+
  theme(
    plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    legend.background = element_rect(fill = "white", color = NA), # White background for legend
    legend.position = "right"
  )

ggsave(paste0(output_fig_ins, "Top_20_Enriched_Pathways.png") , plot = pathways, width = 10, height = 6, dpi = 300)
################# Dot plot for NES and pathways

# Prepare the data
topPathways2 <- fgseaRes %>% 
  arrange(padj) %>% 
  head(n = 20) %>% 
  mutate(pathway = reorder(pathway, NES))

# Create the dot plot
ggplot(topPathways2, aes(x = NES, y = pathway)) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  xlab("Normalized Enrichment Score (NES)") +
  ylab("Pathway") +
  ggtitle("Top 20 Enriched Pathways by NES")


#--------------Inspecting genes in pathways of interest-----------
topPathways2
GOBP_POSITIVE_REGULATION_OF_INTEGRIN_MEDIATED_SIGNALING_PATHWAY
# immunoglobulin_complex_genes <- go[["GOBP_SENSORY_PERCEPTION_OF_CHEMICAL_STIMULUS"]]
# immunoglobulin_complex_genes <- go[["GOBP_POSITIVE_REGULATION_OF_INTEGRIN_MEDIATED_SIGNALING_PATHWAY"]]
# immunoglobulin_complex_genes <- go[["GOBP_GLYCINE_TRANSPORT"]]
immunoglobulin_complex_genes <- go[["GOBP_POSITIVE_REGULATION_OF_INTEGRIN_MEDIATED_SIGNALING_PATHWAY"]]



# Filter the graph_test results to include only genes in this pathway
traj_results_filtered <- traj_results %>%
  filter(gene_short_name %in% immunoglobulin_complex_genes)

# Sort the filtered results by Moran's I statistic to identify the most influential genes
sorted_traj_results <- traj_results_filtered %>%
  arrange(desc(morans_I)) %>%
  dplyr::select(gene_short_name, morans_I, p_value)

# sorted_traj_results <- traj_results_filtered %>%
#   arrange() %>%
#   dplyr::select(gene_short_name, morans_I, p_value)

print(sorted_traj_results)

######################### 

# Visualizing - GOCC_IMMUNOGLOBULIN_COMPLEX genes
features1 = rownames(sorted_traj_results)[1:4]
features2 = rownames(sorted_traj_results)[5:8]
features3 = rownames(sorted_traj_results)[9:12]
FeaturePlot(ins_cells, features = features1, pt.size=1.5)
FeaturePlot(ins_cells, features = features2, pt.size=1.5)
FeaturePlot(ins_cells, features = features3, pt.size=1.5)


FeaturePlot(ins_cells, features = "PCSK1N", pt.size=1.5)
FeaturePlot(ins_cells, features = "ACTB", pt.size=1.5)
FeaturePlot(ins_cells, features = "HSPB1", pt.size=1.5)
FeaturePlot(ins_cells, features = "HLA-H", pt.size=1.5)MTCO3P22
FeaturePlot(ins_cells, features = "MTCO3P22", pt.size=1.5)
FeaturePlot(ins_cells, features = c('BACH2','ALDH1A3'), pt.size=1.5) # Regulators of beta cell dedifferentiation
RBP4, FFAR4/GPR120, ID1, ID2, and ID3



FeaturePlot(ins_cells, features = c("PDX1", "NKX"), pt.size=1.5) #Beta cell markers
FeaturePlot(ins_cells, features = c('BACH2','ALDH1A3'), pt.size=1.5)
FeaturePlot(ins_cells, features = c("NKX6-1", "MYCL"), pt.size=1.5) #Beta cell markers
FeaturePlot(ins_cells, features = c("MYC", "MYCL", "MYCN"), pt.size=1.5) #progenitor cell markings
FeaturePlot(ins_cells, features = "MAFA", pt.size=1.5) #key transcription factor

#-----------------------------From bulk rna-seq------------------------------#
FeaturePlot(ins_cells, features = c("AMY2A", "INS"), pt.size=1.5) #key transcription factor
FeaturePlot(integrated_data, features = c("AMY2A"), pt.size=1.5)


#----------------------integrated gene explortion--------------------------
DimPlot(integrated_data, reduction = "umap", label = TRUE)
FeaturePlot(integrated_data, features = c("PDX1", "NKX6-1"), pt.size=1.5) #Beta cell markers
FeaturePlot(integrated_data, features = c("Ngn3", "Nanog", "OCT4", "L-MYC"))
FeaturePlot(integrated_data, features = c("ALDH1A3"))
FeaturePlot(integrated_data, features = c("MAFA"), pt.size=1.5)


grep("aldh1a3", tolower(rownames(ins_cells)), value = TRUE)
grep("maf", tolower(rownames(ins_cells)), value = TRUE)

#----------------------Searching for branches--------------------------------
unique(integrated_data@meta.data$patient_name)
length(unique(integrated_data@meta.data$patient_name))
