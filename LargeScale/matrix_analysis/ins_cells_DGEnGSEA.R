library(Seurat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(org.Hs.eg.db)
library(fgsea)

integrated_data <- readRDS("C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/intr_labeled.rds")
ins_cells <- readRDS(ins_cells, "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/insulin_cells.rds")
beta_cells <- subset(integrated_data, idents = c("Beta_T2D", "Beta_Healthy"))
output_fig_ins <- "C://Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Figures/SecondSet/ins_cells/ins_vs_beta/"

#--------------------------Normalizing and scaling Data

# Normalize the data
beta_cells <- NormalizeData(beta_cells)
ins_cells <- NormalizeData(ins_cells)

# Find variable features
beta_cells <- FindVariableFeatures(beta_cells)
ins_cells <- FindVariableFeatures(ins_cells)

# Scale the data
beta_cells <- ScaleData(beta_cells)
ins_cells <- ScaleData(ins_cells)

ins_cells@meta.data$condition
# Perform DGEA
dge_results_beta <- FindMarkers(beta_cells, ident.1 = 'Beta_T2D', 
                                ident.2 = 'Beta_Healthy')

ins_cells$ident <- ifelse(ins_cells@meta.data$condition == "Healthy", "Beta_Healthy", "Beta_T2D")
# Set the active identity class to this new column
Idents(ins_cells) <- ins_cells$ident
dge_results_ins <- FindMarkers(ins_cells, ident.1 = "Beta_T2D", 
                               ident.2 = "Beta_Healthy", logfc.threshold = 0.25)

dge_results_beta <- dge_results_beta %>%
  arrange(p_val_adj)
dge_results_ins <- dge_results_ins %>%
  arrange(p_val_adj)


head(dge_results_beta, n=10)
head(dge_results_ins, n=10)

#------------------------------Volcano plot------------------------------------
# Function to create a volcano plot with top 20 gene labels
create_volcano_plot <- function(dge_results, title) {
  # Identify the top 20 genes by adjusted p-value
  top_genes <- dge_results %>%
    arrange(p_val_adj) %>%
    head(20) %>%
    rownames()
  
  # Create the volcano plot
  ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(alpha = 0.4) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-Value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
    geom_text_repel(data = dge_results %>% rownames_to_column("gene") %>% filter(gene %in% top_genes),
                    aes(label = gene),
                    size = 3,
                    max.overlaps = 20)+
    theme(
      plot.background = element_rect(fill = "white", color = NA), # Set plot background to white
      panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
      legend.background = element_rect(fill = "white", color = NA), # White background for legend
      legend.position = "right"
    )
}


# Create volcano plots for both datasets
volcano_beta <- create_volcano_plot(dge_results_beta, "Beta Cells DGEA")
volcano_ins <- create_volcano_plot(dge_results_ins, "Ins Cells DGEA")

# Arrange plots side by side
combined_plot <- plot_grid(volcano_beta, volcano_ins, labels = c("A", "B"))

# Save the combined plot
ggsave(paste0(output_fig_ins, "DGEA_combined_plot.png"), combined_plot, width = 12, height = 6)

#-------------------Identifying the top 20 genes--------------------------------


#-------------------------------GSEA---------------------------------------
# Prepare gene list for GSEA
dge_results_beta <- dge_results_beta %>%
  rownames_to_column("gene") %>%
  mutate(logFC = avg_log2FC) %>%
  arrange(desc(logFC)) %>%
  dplyr::select(gene, logFC)

# Create a named vector
dge_results_beta <- setNames(dge_results_beta$logFC, dge_results_beta$gene)


# Prepare gene list for GSEA
gene_list_ins <- dge_results_ins %>%
  rownames_to_column("gene") %>%
  mutate(logFC = avg_log2FC) %>%
  arrange(desc(logFC)) %>%
  dplyr::select(gene, logFC)

# Create a named vector
gene_list_ins <- setNames(gene_list_ins$logFC, gene_list_ins$gene)


####################GSEA

go <- gmtPathways("C:/Users/jonan/Documents/Tyseq/Data/c5.go.v2023.1.Hs.symbols.gmt")

# Run GSEA
set.seed(06212022)  # Setting a seed for reproducibility
fgseaRes_beta <- fgsea(pathways = go,
                  stats = dge_results_beta,
                  minSize = 6,
                  maxSize = 500,
                  nproc = 1)

fgseaRes_ins <- fgsea(pathways = go,
                  stats = gene_list_ins,
                  minSize = 6,
                  maxSize = 500,
                  nproc = 1)


head(fgseaRes_beta)
head(fgseaRes_ins)


#------------------------------Visualizing up and down--------------------------

choosing_topdown <- function(data){
  # Select top 10 up and down-regulated pathways
  topup <- data %>% filter(ES > 0)
  topup <- topup[order(topup$padj),]
  topup <- topup[1:10,]
  
  topdown <- data %>% filter(ES < 0)
  topdown <- topdown[order(topdown$padj),]
  topdown <- topdown[1:10,]
  
  top_go <- rbind(topup, rev(topdown))
  
  # Clean up pathway names
  top_go$pathway <- top_go$pathway %>% str_replace("KEGG+?_", "") %>% str_replace_all("_", " ")
  
  return(top_go)
}

top_beta <- choosing_topdown(fgseaRes_beta)
top_ins <- choosing_topdown(fgseaRes_ins)


# Function to create and return the plot
# Function to create and return the plot with a legend
create_pathway_figure <- function(data, title){
  # Choose top pathways
  top <- choosing_topdown(data)
  
  # Summarize pathway information and create negative log p-value variable for graphing
  pathg <- top %>% filter(pval <= 0.1) 
  pathg <- pathg %>% mutate(neglogpvalue = -log10(pval))
  
  if (nrow(pathg) >= 1) {
    # Graph pathways by p-value
    pathfig <- ggplot(pathg, aes(x = reorder(pathway, neglogpvalue), y = neglogpvalue)) +
      geom_bar(aes(fill = factor(ifelse(ES < 0, "Downregulated", "Upregulated"))), stat = "identity") +
      coord_flip() +
      scale_x_discrete(name = "Pathways Associated with Diabetes") +
      scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "black"), name = "Regulation") + 
      ylab("-log(p value)") +
      ggtitle(title) +
      theme(axis.text.x = element_text(face = "bold", size = 10, angle = 0),
            axis.text.y = element_text(face = "bold", size = 10, angle = 0),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 10))
    
    return(pathfig)
  } else {
    return(NULL)
  }
}

top_beta_plot <- create_pathway_figure(top_beta, "Top pathways")
top_ins_plot <- create_pathway_figure(top_ins, "Top pathways")

ggsave(paste0(output_fig_ins, "top beta.png"), top_beta_plot,  width = 12, height = 8, dpi = 300)
ggsave(paste0(output_fig_ins, "top ins.png"), top_ins_plot,  width = 12, height = 8, dpi = 300)











