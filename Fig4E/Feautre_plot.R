library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)
library(Seurat)
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


options(repr.plot.width=12, repr.plot.height=10)

# Create each subplot and add labels
p1 <- FeaturePlot(bm, features = "Ly6g", reduction = "tsne") + 
  # Use color gradient similar to example
  scale_color_gradientn(
    colors = c("lightgrey", "#f37907", "#f40404", "#000000"),
    breaks = c(100, 200, 300, 400),
    name = "No. of\nneighbours"
  )+   theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5)
  )

p2 <- FeaturePlot(bm, features = "Csf1r", reduction = "tsne") + 
  # Use color gradient similar to example
  scale_color_gradientn(
    colors = c("lightgrey", "#f37907", "#f40404", "#000000"),
    breaks = c(100, 200, 300, 400),
    name = "No. of\nneighbours"
  )+   theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5)
  )

p3 <- FeaturePlot(bm, features = "Arg2", reduction = "tsne") + 
  # Use color gradient similar to example
  scale_color_gradientn(
    colors = c("lightgrey", "#f37907", "#f40404", "#000000"),
    breaks = c(100, 200, 300, 400),
    name = "No. of\nneighbours"
  )+   theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5)
  )

p4 <- FeaturePlot(bm, features = "Il1b", reduction = "tsne") + 
  # Use color gradient similar to example
  scale_color_gradientn(
    colors = c("lightgrey", "#f37907", "#f40404", "#000000"),
    breaks = c(100, 200, 300, 400),
    name = "No. of\nneighbours"
  )+   theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5)
  )

# Combine subplots
combined_plot <- (p1 | p2) / (p3 | p4)



pdf("./Feature.pdf", width = 12, height = 12)
print(combined_plot)
dev.off()