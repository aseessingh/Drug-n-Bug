library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(tidyr)
# Read the CSV file (make sure to provide the correct path)
effect_size_data <- read_csv("Chemical Annotations Integrated modD greater than 0.2.csv", show_col_types = FALSE)

# Clean the data: Remove rows with missing values and ensure values are within the scale range
effect_size_data <- effect_size_data %>% select(`Spearman Rho serum metformin`,`Spearman Rho HbA1c`,CHEMICAL_NAME,Source)
effect_size_data <- effect_size_data %>% mutate_all(~replace_na(., "Unknown"))
# Create the scatter plot
plot <- ggplot(effect_size_data, aes(x = `Spearman Rho serum metformin`, 
                                    y = `Spearman Rho HbA1c`, 
                                    color = Source)) +
  geom_point(size = 3, alpha = 0.6) +  # Use transparency
  geom_text_repel(aes(label = CHEMICAL_NAME), 
                  size = 2.8, 
                  vjust = -1, 
                  hjust = 1, 
                  max.overlaps = Inf,  # Allow unlimited overlaps
                  nudge_x = 0.05,      # Nudge labels slightly to the right
                  nudge_y = 0.05) +     # Nudge labels slightly up
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown", 
                              "pink", "yellow", "cyan", "magenta", "gray", "black", 
                              "darkblue", "darkgreen", "lightblue", "lightgreen", "darkred", "darkorange")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Horizontal axis at y = 0
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(title = "Scatter Plot of Common Associated Metabolites",
       x = "Spearman Rho of Serum Metformin",
       y = "Spearman Rho of HbA1c") +
  theme_minimal(base_size = 11) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  # Set the limits for the x and y axes
  scale_x_continuous(limits = c(-1, 1)) +  # Adjust these limits as needed
  scale_y_continuous(limits = c(-1, 1))
# Save the plot as a PNG file
ggsave("plot.png", 
       plot = plot, width = 10, height = 8, dpi = 600)
