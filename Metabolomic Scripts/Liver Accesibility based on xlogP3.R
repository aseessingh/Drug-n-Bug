library(ggplot2)
library(ggrepel)
library(readr)

# Read the CSV file (make sure to provide the correct path)
effect_size_data <- read_csv("adj_annotated_logp.csv", show_col_types = FALSE)
effect_size_data$liver_access <- ifelse(effect_size_data$xlogP3 > 2.5, 
                                        "can enter liver", 
                                        "cannot enter liver")

plot <- ggplot(effect_size_data, aes(x = `Effect Size of Serum Metformin`, 
                                     y = `Effect Size of HbA1c`, 
                                     color = liver_access, 
                                     shape = liver_access)) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_text_repel(aes(label = CHEMICAL_NAME),
                  size = 2.8,
                  max.overlaps = 1000,
                  force = 5,
                  force_pull = 1,
                  box.padding = 0.8,
                  point.padding = 0.4,
                  segment.size = 0.2,
                  segment.alpha = 0.4,
                  min.segment.length = 0.1,
                  nudge_y = 0.01) +
  
  scale_color_manual(values = c("can enter liver" = "#1b9e77", 
                                "cannot enter liver" = "#d95f02")) +
  scale_shape_manual(values = c("can enter liver" = 16, 
                                "cannot enter liver" = 17)) +
  
  labs(title = "Effect Size Plot by Liver Accessibility (xlogP3)",
       x = "Effect Size of Serum Metformin",
       y = "Effect Size of HbA1c",
       color = "Liver Accessibility",
       shape = "Liver Accessibility") +
  
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
ggsave("Adjusted_Association_Common_Plot.png", plot = plot, width = 8, height = 7, dpi = 600)
