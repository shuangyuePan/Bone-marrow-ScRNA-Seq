library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# Define promoter region (chr14: 67,619,000-67,620,998)
promoter_region <- GRanges(
  seqnames = "14", 
  ranges = IRanges(start = 67622739, end = 67623450)
)

# Input files (relative paths from FigS4I directory to Fig3E where BW files are)
# Assuming BW files are in Fig3E as per previous steps
ko_bw_path <- "../Fig3E/GSE237222_KO_merged_ratio.bw"
nc_bw_path <- "../Fig3E/GSE237222_NC_merged_ratio.bw"

# Import data
ko_coverage <- import(ko_bw_path, which = promoter_region)
nc_coverage <- import(nc_bw_path, which = promoter_region)

# Add metadata
ko_coverage$sample <- "KO"
nc_coverage$sample <- "NC"

# Combine data
combined_data_df <- as.data.frame(c(ko_coverage, nc_coverage))

# Define target region and sequence
target_xmin <- 67623165
target_xmax <- 67623169
target_seq  <- "ATTAAA" 

# Text box position
box_x <- (target_xmin + target_xmax) / 2 
box_y <- max(combined_data_df$score) * 0.15

# Highlight region
highlight_region <- data.frame(
  xmin = 67623165, 
  xmax = 67623169, 
  ymin = -Inf, 
  ymax = Inf
)

# Plot
p <- ggplot(combined_data_df, aes(x = start, y = score, color = sample)) +
  # 1. Highlight background
  geom_rect(
    data = highlight_region,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "orange", alpha = 0.15, inherit.aes = FALSE
  ) +
  # 2. Plot smooth lines
  geom_line(linewidth = 1.5, alpha = 0.85) +
  
  # 3. Add text annotation
  annotate( 
    "text", x = 67623100, y = 28, # Adjusted Y position
    label = "4-fold enrichment\nin KO vs NC", color = "darkred", size = 4, fontface = "bold" 
  ) + 
  # 4. Add double-headed arrow
  annotate( 
    "segment", x = 67623160, xend = 67623160, 
    y = 20, yend = 45, # Adjusted based on visual estimation
    arrow = arrow(type = "closed", length = unit(0.2, "cm"), ends = "both"), 
    color = "darkred", linewidth = 1.2 
  ) +
  
  # 5. Add text box and leads
  # Left lead
   annotate("segment", 
            x = box_x - 10,  
            y = box_y,       
            xend = target_xmin, 
            yend = 0,           
            color = "darkred", linewidth = 0.8) +
   # Right lead
   annotate("segment", 
            x = box_x + 10,     
            y = box_y, 
            xend = target_xmax, 
            yend = 0, 
            color = "darkred", linewidth = 0.8) +
  # Text box
  annotate("label", 
           x = box_x, 
           y = box_y, 
           label = target_seq, 
           fill = "white",      
           color = "darkred",     
           size = 4, 
           fontface = "bold",
           label.size = 0.5) +

  # 6. Colors and Styling
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
    legend.position = c(0.15, 0.9), # Move legend to top-left inside plot
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(), # Hide x-axis text
    axis.ticks.x = element_line(),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15))
  ) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_line())

# Save output
output_file <- "FigS4F_ChIP_seq_coverage.pdf"
ggsave(output_file, plot = p, width = 6, height = 5) # Adjusted size to be more square-like
print(paste("Plot saved to:", output_file))
