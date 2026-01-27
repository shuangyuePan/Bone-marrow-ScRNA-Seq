library(DESeq2)
library(data.table)
library(ggplot2)

input_file <- "./counts_anno.xls"
out_pdf_pca <- "FigS4A_PCA.pdf"

dt <- data.table::fread(input_file, sep = "\t", header = TRUE)
all_cols <- names(dt)
sample_cols <- all_cols[grep("KO|NC", all_cols, ignore.case = TRUE)]

count_matrix <- as.matrix(dt[, ..sample_cols])
rownames(count_matrix) <- dt$id
count_matrix <- round(count_matrix)
count_matrix <- count_matrix[!is.na(rowMeans(count_matrix)) & is.finite(rowMeans(count_matrix)), ]
count_matrix[count_matrix < 0] <- 0

sample_info <- data.frame(
  group = ifelse(grepl("KO", sample_cols, ignore.case = TRUE), "KO", "NC"),
  row.names = sample_cols
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ group
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

pcaData_vst <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

color_scheme <- c("KO" = rgb(234/255, 156/255, 157/255),
                  "NC" = rgb(186/255, 213/255, 231/255))

get_axis_limits <- function(pca_data) {
  x_range <- range(pca_data$PC1)
  y_range <- range(pca_data$PC2)
  max_range <- max(abs(c(x_range, y_range)))
  return(c(-max_range * 1.1, max_range * 1.1))
}

axis_limits_vst <- get_axis_limits(pcaData_vst)

pdf(out_pdf_pca, width = 6, height = 6)

p1 <- plotPCA(vsd, intgroup = "group") + 
  ggtitle("PCA of samples - VST transformation") +
  theme_minimal() +
  scale_color_manual(values = color_scheme) +
  xlim(axis_limits_vst) +
  ylim(axis_limits_vst) +
  coord_fixed(ratio = 1)
print(p1)

dev.off()
