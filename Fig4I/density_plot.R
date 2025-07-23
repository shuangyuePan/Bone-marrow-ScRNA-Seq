library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)
library(Matrix)
library(parallel)
library(tidyverse)
library(patchwork)
library(harmony)
library(ggpubr)
library(magrittr)
library(FNN) # For distance-based neighbor search
library(ggforce)
library(wesanderson)
library(qs)

bm <- qread("../bm.qs")

#-----------------------------------------------------------
# 1. Extract data
#-----------------------------------------------------------


tsne_coords <- Embeddings(bm, reduction = "tsne")

# Create base dataframe
plot_data <- data.frame(
  tSNE_1 = tsne_coords[,1],
  tSNE_2 = tsne_coords[,2],
  celltype = bm$celltype_final, # Use cell type information
  group = bm$group          # Use group information (KO/NC)
)



#-----------------------------------------------------------
# 2. Calculate neighbor counts by group
#-----------------------------------------------------------


# Calculate neighbor counts separately for each group
unique_groups <- unique(plot_data$group)
neighbor_counts <- numeric(nrow(plot_data))

# Calculate neighbor counts by group
for (current_group in unique_groups) {
  # Get indices for current group
  group_idx <- which(plot_data$group == current_group)
  
  # Extract tSNE coordinates for current group
  group_coords <- tsne_coords[group_idx, ]
  
  # Use FNN package to calculate neighbors
  k_max <- min(600, length(group_idx) - 1)  # Ensure k is not larger than group cell count-1
  # Set search radius
  search_radius <- 2.66
    
  if (k_max > 0) {
    nn_result <- get.knnx(data = group_coords, query = group_coords, k = k_max)
    
    # Calculate number of neighbors within radius
    group_neighbor_count <- rowSums(nn_result$nn.dist <= search_radius) - 1
    
    # Store results in main dataframe
    neighbor_counts[group_idx] <- group_neighbor_count
  }
}

# Add neighbor counts to dataframe
plot_data$neighbors <- neighbor_counts
bm$neighbors <- neighbor_counts # Also save to Seurat object




#-----------------------------------------------------------
# 3. Prepare data for plotting
# Define cell types to highlight
 highlight_celltypes <- c("PMN-MDSC")

# Create mask for highlighted cells
highlighted_cells <- plot_data$celltype %in% highlight_celltypes
plot_data$highlight <- highlighted_cells

# Extract data for highlighted cells
highlight_data <- plot_data[highlighted_cells, ]

# Create gray background data (all cells)
background_data <- plot_data
background_data$neighbors <- 0  # Set to 0 to display as gray


#-----------------------------------------------------------
# 4. Statistical testing and result annotation
#-----------------------------------------------------------
# Create contingency table
contingency_table <- table(
  Group = plot_data$group,
  Highlight = plot_data$celltype %in% highlight_celltypes
)

# Fisher's exact test
cat("\nFisher's exact test results:\n")
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
print(fisher_result)

# Format p-value
p_value <- ifelse(fisher_result$p.value < 0.001, "p_val < 0.001", 
                 sprintf("p_val = %.3f", fisher_result$p.value))

# Create annotation data
annotation_data <- data.frame(
  label = p_value,
  x = mean(range(plot_data$tSNE_1)),  # Center on X axis
  y = max(plot_data$tSNE_2) * 1.1     # Above the top of Y axis
)

#-----------------------------------------------------------
# 5. Create group comparison plot
#-----------------------------------------------------------
# Create facet plot
p_facet <- ggplot() +
  # Add gray background points
  geom_point(data = background_data, 
             aes(x = tSNE_1, y = tSNE_2), 
             color = "lightgrey", 
             size = 1, 
             alpha = 0.5) +
  # Add highlighted cells, colored by neighbor count
  geom_point(data = highlight_data, 
             aes(x = tSNE_1, y = tSNE_2, color = neighbors), 
             size = 1)  +
  # Use color gradient similar to example
  scale_color_gradientn(
    colors = c("#FFFFCC", "#f37907", "#f40404", "#000000"),
    limits = c(0, 380),
    breaks = c(100, 200, 300),
    name = "No. of\nneighbours"
  ) +
  # Facet by group
  facet_wrap(~ group, scales = "fixed") +
  # Add global title
  ggtitle(paste( p_value) ) +
  # Expand plot area and disable clipping
  coord_cartesian(ylim = c(min(plot_data$tSNE_2), max(plot_data$tSNE_2) * 1.2),
                  clip = "off") +
  # Beautify plot
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1,
                              ),  # t=top margin, b=bottom margin
    panel.grid = element_blank(),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20,face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10,face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  labs(x = "t-SNE 1", y = "t-SNE 2")



pdf("./PMN.All.tsne.pdf", width = 12, height = 6.5)
print(p_facet)
dev.off()