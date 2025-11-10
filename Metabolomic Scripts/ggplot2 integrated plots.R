# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(readr)

# Read the CSV file (make sure to provide the correct path)
effect_size_data <- read_csv("adj_annotated_common_metabolites.csv", show_col_types = FALSE)

# Create the scatter plot
plot <- ggplot(effect_size_data, aes(x = `Effect Size of Serum Metformin`, 
                                     y = `Effect Size of HbA1c`, 
                                     color = Source, 
                                     shape = Source)) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_text_repel(aes(label = CHEMICAL_NAME),
                  size = 2.8,                         # Smaller text
                  max.overlaps = 1000,                # Allow many overlaps to be checked
                  force = 5,                          # Stronger repulsion
                  force_pull = 1,                     # Prevent labels from snapping back
                  box.padding = 0.8,                  # More space around label
                  point.padding = 0.4,
                  segment.size = 0.2,
                  segment.alpha = 0.4,
                  min.segment.length = 0.1,
                  nudge_y = 0.01) +                   # Slight nudge can reduce stacking
  
  scale_color_brewer(palette = "Dark2") +             # More readable palette
  labs(title = "Scatter Plot of Common Associated Metabolites",
       x = "Effect Size of Serum Metformin",
       y = "Effect Size of HbA1c") +
  theme_minimal(base_size = 11) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Save the plot as a PNG file
ggsave("Adjusted_Association_Common_Plot.png", plot = plot, width = 8, height = 7, dpi = 600)
