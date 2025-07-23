library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# Define promoter region (using chromosome "14", coordinate range 67,619,000-67,620,998)
promoter_region <- GRanges(
  seqnames = "14", 
  ranges = IRanges(start = 67619000, end = 67621000)
)

# Import bw files and add sample labels
ko_coverage <- import("../Fig3E/GSE237222_KO_merged_ratio.bw", which = promoter_region)
nc_coverage <- import("../Fig3E/GSE237222_NC_merged_ratio.bw", which = promoter_region)
ko_coverage$sample <- "KO"
nc_coverage$sample <- "NC"

# Combine data
combined_data <- c(ko_coverage, nc_coverage)

# Convert GenomicRanges object to dataframe
combined_data_df <- as.data.frame(combined_data)

# Define KO significantly high expression region (assumed in the 67,619,500-67,620,000 range)
highlight_region <- data.frame(
  xmin = 67619800, 
  xmax = 67620200, 
  ymin = -Inf, 
  ymax = Inf
)

# Plot coverage and add annotations
ggplot(combined_data_df, aes(x = start, y = score, color = sample)) +
  # 1. Highlight differential region background
  geom_rect(
    data = highlight_region,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "orange", alpha = 0.15, inherit.aes = FALSE
  ) +
  # 2. Draw smoothed lines (optimize visual smoothness with line width and transparency)
  geom_line(size = 1.5, alpha = 0.85) +
  # 3. Add fold change annotation
  annotate(
    "text", x = 67619350, y = max(combined_data_df$score) * 0.5,
    label = "10-fold enrichment\nin KO vs NC", color = "darkred", size = 4, fontface = "bold"
  ) +
  # 4. Add pointing arrow
  annotate(
    "segment", x = 67619850, xend = 67619850,
    y = max(subset(combined_data_df, sample == "KO")$score) * 0.95,
    yend = max(subset(combined_data_df, sample == "NC")$score) * 1.2,
    arrow = arrow(type = "closed", length = unit(0.2, "cm"), ends = "both"),
    color = "darkred", size = 1.2
  ) +
  # 5. Adjust color theme and axes
  scale_color_manual(
    values = c("KO" = "#BAD5E7", "NC" = "#EA9C9D"),
    labels = c("KO" = expression(paste("ARID1A"^"KO")), "NC" = expression(paste("ARID1A"^"NC"))),
    limits = c("NC", "KO")  # Set legend order, NC on top
  ) +
  labs(
    title = "ChIP-seq Coverage at ARG2 Promoter",
    x = "Relative Position", 
    y = "ChIP-seq Signal",
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(r = 10)),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(face = "italic", margin = margin(t = 15)),
    axis.title.y = element_text(face = "italic", margin = margin(r = 15))
  )

# Save plot as PDF

ggsave("Arg2_peak_plot.pdf")