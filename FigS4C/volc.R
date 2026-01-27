library(ggplot2)

input_file <- "./KO-vs-NC-all.gene.xls"
output_file <- "FigS4C_KO-vs-NC_volcano_plot.pdf"

ko_color <- rgb(234/255, 156/255, 157/255)
nc_color <- rgb(186/255, 213/255, 231/255)

data <- read.table(input_file, header = TRUE, sep = "\t", quote = "")
colnames(data)[colnames(data) == "gene_id"] <- "gene"
if("q-value" %in% colnames(data)) {
  colnames(data)[colnames(data) == "q-value"] <- "q.value"
}

data$color <- "grey"
data$color[data$log2FoldChange > 0.5 & data$q.value < 0.05] <- "red"
data$color[data$log2FoldChange < -0.5 & data$q.value < 0.05] <- "blue"

x_max <- max(abs(data$log2FoldChange), na.rm = TRUE) * 1.2
x_lim <- c(-x_max, x_max)

p <- ggplot(data, aes(x = log2FoldChange, y = -log10(q.value), color = color)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("red" = ko_color, "blue" = nc_color, "grey" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "log2(Fold Change)", y = "-log10(q-value)", title = "Volcano Plot (KO vs NC)") +
  xlim(x_lim)

pdf(output_file, width = 6, height = 6)
print(p)
dev.off()
