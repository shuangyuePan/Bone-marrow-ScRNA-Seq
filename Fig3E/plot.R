library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# Define promoter region (chr14: 67,619,000-67,622,000)
promoter_region <- GRanges(
  seqnames = "14", 
  ranges = IRanges(start = 67619000, end = 67622000)
)

# Input files (relative paths to current directory)
ko_bw_path <- "./GSE237222_KO_merged_ratio.bw"
nc_bw_path <- "./GSE237222_NC_merged_ratio.bw"

# Import data
ko_coverage <- import(ko_bw_path, which = promoter_region)
nc_coverage <- import(nc_bw_path, which = promoter_region)

# Add metadata
ko_coverage$sample <- "KO"
nc_coverage$sample <- "NC"

# Combine data
combined_data_df <- as.data.frame(c(ko_coverage, nc_coverage))

# Define highlight region for 10-fold enrichment
highlight_region <- data.frame(
  xmin = 67619800, 
  xmax = 67620200, 
  ymin = -Inf, 
  ymax = Inf
)

# Plot
p <- ggplot(combined_data_df, aes(x = start, y = score, color = sample)) +
  
  # 1. Highlight background region
  geom_rect(
    data = highlight_region,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "orange", alpha = 0.15, inherit.aes = FALSE
  ) +
  
  # 2. Plot smooth lines
  geom_line(linewidth = 1.5, alpha = 0.85) +
  
  # 3. Add text annotation
  annotate(
    "text", x = 67619350, y = 25, # Adjusted y position based on plot scale
    label = "10-fold enrichment in\nARID1AKO vs ARID1ANC", 
    color = "darkred", size = 4, fontface = "bold", hjust = 0
  ) +
  
  # 4. Add double-headed arrow
  annotate(
    "segment", x = 67619900, xend = 67619900,
    y = 8, yend = 48, # Approximate positions based on visual inspection
    arrow = arrow(type = "closed", length = unit(0.2, "cm"), ends = "both"),
    color = "darkred", linewidth = 1.2
  ) +
  
  # 5. Styling and Labels
  scale_color_manual(
    values = c("KO" = "#BAD5E7", "NC" = "#EA9C9D"),
    labels = c("KO" = expression(paste("ARID1A"^"KO")), "NC" = expression(paste("ARID1A"^"NC"))),
    limits = c("NC", "KO")
  ) +
  labs(
    x = "Relative Position", 
    y = "ChIP-seq Signal",
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.6), # Inside plot legend
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(), # Hide x-axis text as per example image style (or keep if needed)
    axis.ticks.x = element_line(),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  # Remove x-axis labels to match "Relative Position" style usually without genomic coords
  theme(axis.text.x = element_blank(), axis.ticks.x = element_line())

# Save output
output_file <- "Fig3E_CHIP_seq_coverage.pdf"
ggsave(output_file, plot = p, width = 6, height = 5) # Adjusted size to be more square-like
print(paste("Plot saved to:", output_file))
