library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(stringr)
library(forcats)
library(qs)

# Input data path (adjust as needed)
input_qs <- "../../0720_paper/bm.qs"
output_pdf <- "FigS10A_Enrichment_Dotplot.pdf"

# Load Seurat object
if (file.exists(input_qs)) {
  seu <- qread(input_qs)
} else {
  stop("Input Seurat object not found at ", input_qs)
}

# Define cell type and comparison
seu$celltype_final <- as.character(seu$celltype_final) # Ensure character vector
seu$cell_sur <- ifelse(seu$celltype_final == "PMN-MDSC", "case", "control")

# Find markers (KO PMN vs ALL or similar logic based on context)
# Assuming 'group' column exists for KO/NC distinction within PMN-MDSC if that was the intent,
# OR finding markers for PMN-MDSC vs others. 
# Based on previous code snippet logic: `ident.1` = "case", `ident.2` = "control"
Idents(seu) <- "cell_sur"
markers <- FindMarkers(seu, ident.1 = "case", ident.2 = "control", max.cells.per.ident = 2000)
markers$gene <- rownames(markers)

# Filter significant markers
filtered <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
filtered_symbols <- unique(filtered$gene)

# Convert Symbols to Entrez IDs
gene_df <- bitr(filtered_symbols, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Mm.eg.db)

# 1. GO Enrichment (Biological Process)
ego <- enrichGO(gene          = gene_df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# Select specific GO terms for FigS10A
target_go_terms <- c("amine metabolic process", "biogenic amine metabolic process")
go_data <- as.data.frame(ego) %>%
  filter(Description %in% target_go_terms) %>%
  mutate(Type = "GO: BP") %>%
  select(Description, GeneRatio, p.adjust, Count, Type)

# 2. Reactome Enrichment
rpa <- enrichPathway(gene = gene_df$ENTREZID,
                     organism = "mouse",
                     pvalueCutoff = 0.5,
                     readable = TRUE)

# Select specific Reactome terms for FigS10A
target_reactome_terms <- c("Metabolism of polyamines")
ra_data <- as.data.frame(rpa) %>%
  filter(Description %in% target_reactome_terms) %>%
  mutate(Type = "Reactome") %>%
  select(Description, GeneRatio, p.adjust, Count, Type)

# Combine results
plot_data <- rbind(go_data, ra_data)

# Calculate GeneRatio as numeric
plot_data$GeneRatioNum <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

# Order descriptions
plot_data$Description <- factor(plot_data$Description, levels = c("Metabolism of polyamines", "Biogenic amine metabolic process", "Amine metabolic process"))
plot_data$Type <- factor(plot_data$Type, levels = c("GO: BP", "Reactome"))

# Plotting
p <- ggplot(plot_data, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradientn(
    colours = c("#FF0000", "#FF00FF", "#0000FF"), 
    limits = c(0, 0.1),
    name = "p.adjust",
    guide = guide_colorbar(order = 1)
  ) +
  scale_size_continuous(range = c(5, 10), name = "Gene Count") +
  facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    title = "Enrichment of Amine & Polyamine Metabolism",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text.y = element_text(angle = 270, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    strip.background = element_rect(fill = "grey95")
  )

# Save output
ggsave(output_pdf, p, width = 8, height = 5)
print(paste("Plot saved to:", output_pdf))
