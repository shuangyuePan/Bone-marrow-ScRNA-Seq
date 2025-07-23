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

prop = bm@meta.data %>% dplyr::group_by(celltype_final, group)  %>%  dplyr::summarise(num = dplyr::n()) 

prop$group <- prop$group %>% factor( levels =c("KO","NC"))

combined_plot <- prop %>%
  ggplot( aes(x = celltype_final, y = num, fill = group))+
  theme_pubr() +
  geom_bar(stat = "identity", position = "fill") +
  geom_hline(yintercept = 0.5, linetype = 2, size = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Relative cell proportion") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1), # Adjust x-axis text
    axis.text.y = element_text(size = 12, face = "bold"), # Adjust y-axis text
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold',margin = margin(r=20)),
    legend.title = element_blank()
  ) 


pdf("./cellRatio.pdf", width = 12, height = 12)
print(combined_plot)
dev.off()