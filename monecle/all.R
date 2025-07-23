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
library(monocle)
library(BiocGenerics)

##########################################################################################
bm <- qread("../bm.qs")
# Select more cell types for a broader analysis
bm.neu.2 = subset(bm, celltype %in% c("Early Myeloid", "Neutrophils", "pre-MDSC", "PMN-MDSC"))

# Set random seed for reproducibility and downsample to 800 cells
set.seed("20240802")
bm.neu.2 = subset(bm.neu.2, downsample = 800)

# Visualize cell types
p_dim <- DimPlot(bm.neu.2, group.by = "celltype")
ggsave("celltype_dimplot.pdf", p_dim, width=8, height=6)

# Get expression matrix, metadata, and gene annotations
expression_matrix = bm.neu.2@assays$RNA$counts
cell_metadata = data.frame(bm.neu.2@meta.data)
cell_metadata$batch = factor(bm.neu.2$orig.ident)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")

# Create Monocle2 cell dataset object
cds_neu_m2 <- newCellDataSet(as(expression_matrix, "sparseMatrix"),
                          phenoData  = AnnotatedDataFrame(cell_metadata),
                          featureData  = AnnotatedDataFrame(gene_annotation),
                          expressionFamily=negbinomial.size()
)

# Estimate size factors and dispersions (technical parameters)
cds_neu_m2 <- estimateSizeFactors(cds_neu_m2)
cds_neu_m2 <- estimateDispersions(cds_neu_m2)

# Check feature data
fData(cds_neu_m2)

# Detect which genes are expressed in cells
cds_neu_m2 <- detectGenes(cds_neu_m2, min_expr = 0.1)
print(head(fData(cds_neu_m2)))

# Repeat gene detection (seems redundant)
cds_neu_m2 <- detectGenes(cds_neu_m2, min_expr = 0.1)
print(head(fData(cds_neu_m2)))

# Keep only genes expressed in at least 40 cells
expressed_genes <- row.names(subset(fData(cds_neu_m2),
                                    num_cells_expressed >= 40))

# Check phenotype data
print(head(pData(cds_neu_m2)))

# Find genes that differ by cell type
diff_test_res <- differentialGeneTest(cds_neu_m2[expressed_genes,],
                                      fullModelFormulaStr = "~celltype",
                                      cores = 30
)

# Find genes that differ between knockout (KO) and wildtype (normal control, NC)
diff_test_res1 <- differentialGeneTest(cds_neu_m2[expressed_genes,],
                                      fullModelFormulaStr = "~group",
                                      cores = 30
)

# Select significant genes from both tests
ordering_genes <- row.names(subset(diff_test_res, qval < 0.001))
ordering_genes1 <- row.names(subset(diff_test_res1, qval < 0.01))

# Get dispersion table for genes
disp_table <- dispersionTable(cds_neu_m2)

# Find intersection of genes that differ by both cell type and treatment group
used_genes = intersect(ordering_genes1, ordering_genes)

# Keep only genes with positive empirical dispersion
used_genes = used_genes[used_genes %in% disp_table[disp_table$dispersion_empirical >0 ,"gene_id"]]

# Set the genes used for ordering cells along trajectory
cds_neu_m2 <- setOrderingFilter(cds_neu_m2, used_genes)

# Plot genes used for ordering
p_ordering <- plot_ordering_genes(cds_neu_m2)
ggsave("ordering_genes_plot.pdf", p_ordering, width=8, height=6)

# Reduce dimensions using DDRTree (specific to Monocle2)
cds_neu_m2 <- reduceDimension(cds_neu_m2, max_components = 2,
                           reduction_method = 'DDRTree',
                           verbose = T)

# Order cells along the trajectory
cds_neu_m2 <- orderCells(cds_neu_m2)

# Load ggpubr for publication-quality plots

# Create a custom scatter plot of the trajectory dimensions
p_scatter <- p1$data |>
  ggplot(aes(x= data_dim_1, y = data_dim_2)) +
  theme_pubr() +
  geom_point(aes(col = celltype))
ggsave("trajectory_scatter.pdf", p_scatter, width=8, height=6)

# Define a function to find the starting state for trajectory (root state)
# The root state is typically where undifferentiated cells are found
GM_state <- function(cds, colnames, naive_names){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)[, "State"], pData(cds)[, colnames])[,naive_names]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

# Look at the distribution of cell types across trajectory states
table(pData(cds_neu_m2)[, "State"], pData(cds_neu_m2)[, "celltype"])

# Set the root state to where most Early Myeloid cells are found
cds_neu_m2 <- orderCells(cds_neu_m2, root_state = GM_state(cds_neu_m2, colnames = "celltype",
                                                     naive_names = "Early Myeloid"))



## Fig4F 

# 保存状态轨迹图
p_state <- plot_cell_trajectory(cds_neu_m2, color_by = "State") +
  scale_color_manual(values = my36colors[5:30]) +
  theme_pubr()
ggsave("trajectory_state.pdf", p_state, width=10, height=5.5)

# 保存细胞类型轨迹图
p_celltype <- plot_cell_trajectory(cds_neu_m2, color_by = "celltype") +
  scale_color_manual(values = c("#F3B1A0", "#D6E7A3", "#E5D2DD", "#F1BB72")) +
  theme_pubr()
ggsave("trajectory_celltype.pdf", p_celltype, width=10, height=5.5)

# 保存分组轨迹图
p_group <- plot_cell_trajectory(cds_neu_m2, color_by = "group") +
  scale_color_manual(values = c("#FC8D62", "#4DBBD5")) +
  theme_pubr()
ggsave("trajectory_group.pdf", p_group, width=10, height=5.5)



##Fig4G

# Create a bar plot showing cell type composition in each state
p1 = bm.neu.2@meta.data |>
  ggplot(aes(x = State, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#F3B1A0", "#D6E7A3", "#E5D2DD", "#F1BB72") )+
  scale_y_continuous(labels = scales::percent) +
  labs(y = "% cells") +
  theme_pubr(legend = "right")

# Plot showing proportion of each treatment group in each state
p_group_comp <- bm.neu.2@meta.data |>
  ggplot(aes(x = State, fill = group)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.5) +
  scale_fill_manual(values = c("#FC8D62", "#4DBBD5")) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "% cells") +
  theme_pubr(legend = "right")
ggsave("group_composition_state.pdf", p_group_comp, width=8, height=6)

# Create a custom bar plot with specific cell counts
p2 = data.frame( State = c(1:5, 1:5),
            group = rep(c("KO", "NC"), each = 5),
            num = c(441, 63, 1211, 123, 66, 515, 48, 400, 245, 88)
            ) |>
  ggplot(aes(x = State, y = num, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.5) +
  scale_fill_manual(values = c("#FC8D62", "#4DBBD5"))+
  labs(y = "% cells") +
  theme_pubr(legend = "right")

# 保存细胞类型组成图
p_celltype_comp <- bm.neu.2@meta.data |>
  ggplot(aes(x = State, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#F3B1A0", "#D6E7A3", "#E5D2DD", "#F1BB72")) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "% cells") +
  theme_pubr(legend = "right")
ggsave("celltype_composition.pdf", p_celltype_comp, width=8, height=6)

# 保存硬编码分组组成图
p_hardcoded_group <- data.frame(State = c(1:5, 1:5),
            group = rep(c("KO", "NC"), each = 5),
            num = c(441, 63, 1211, 123, 66, 515, 48, 400, 245, 88)) |>
  ggplot(aes(x = State, y = num, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.5) +
  scale_fill_manual(values = c("#FC8D62", "#4DBBD5")) +
  labs(y = "% cells") +
  theme_pubr(legend = "right")
ggsave("hardcoded_group_composition.pdf", p_hardcoded_group, width=8, height=6)


#Fig S9A


  gene_list = c('Arg2','Sat1','Picalm','Malat1','Wfdc17','Ifitm1','Srgn','Cebpb','Cstdc4','Dusp1','Stfa2l1',
              'Cd177','Adpgk','Anxa1','Syne1','Lcn2','Aldh2','Serpinb1a','Camp','Ak2','Lyz2','Arhgdib',
              'Ltf','Golim4','Ngp','Cybb','Anxa3','Abca13','Ifitm6','Lgals3','Ly6c2','Capg','Lims1','S100a8','Pglyrp1','Ly6g','Wfdc21')

p <- plot_genes_branched_heatmap(cds_neu_m2[gene_list,],
                             branch_point = 1,
                            branch_labels = c("Cell state3 ", "Cell state4"),
                             num_clusters = 3,
                             cores = 4,
                             use_gene_short_name = T,
                             show_rownames = T,
                           return_heatmap = TRUE
                           )

ggsave("branched_heatmap.pdf", p$ph_res, width=6, height=8)




###  FigS9E 



plot_genes_branched_pseudotime(cds_neu_m2[c("Arg2", "Sat1"),],
                                branch_point = 1,
                                color_by = "State",
                                ncol = 1)+   scale_color_manual(values = my36colors[5:30])

p_pseudotime <- plot_genes_branched_pseudotime(cds_neu_m2[c("Arg2", "Sat1"),],
                                              branch_point = 1,
                                              color_by = "State",
                                              ncol = 1) +
  scale_color_manual(values = my36colors[5:30])
ggsave("branched_pseudotime.pdf", p_pseudotime, width=10, height=8)


                                