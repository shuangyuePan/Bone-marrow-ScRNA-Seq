library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)

# Update input path to be relative to Fig3B directory
input_file <- "../FigS4C/KO-vs-NC-all.gene.xls"
output_csv <- "Fig3B_metabolism_KEGG_enrichment_top10.csv"
output_pdf <- "Fig3B_metabolism_KEGG_dotplot.pdf"

fc <- 0.5
set.seed(123)

data <- data.table::fread(input_file, sep = "\t", header = TRUE, fill = TRUE, data.table = FALSE)

filtered_data <- data %>% filter(abs(log2FoldChange) > fc & `q-value` < 0.05)
diff_genes <- filtered_data$gene_id

gene_df <- bitr(diff_genes, 
                fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)

kegg_result <- enrichKEGG(
    gene = gene_df$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.2,
    qvalueCutoff = 0.3
)

kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

metabolism_result <- as.data.frame(kegg_result) %>%
    filter(category == "Metabolism") %>%
    arrange(p.adjust) %>%
    head(10)

write.csv(metabolism_result, output_csv, row.names = FALSE)

metabolism_result$Description <- factor(metabolism_result$Description, 
                                     levels = rev(metabolism_result$Description))

pdf(output_pdf, width = 12, height = 8)
print(ggplot(metabolism_result, 
            aes(x = -log10(p.adjust), y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    theme_minimal() +
    theme(text = element_text(family = "sans")) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "-log10(p.adjust)", y = "Pathway",
         title = sprintf("Top 10 Metabolism Pathways (|log2FC| > %.1f)", fc)))
dev.off()
