library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggnewscale)
library(patchwork)
results <- as.data.frame(read_csv("top 10 loadings with annotation.csv",show_col_types = FALSE))
drug_data <- as.data.frame(read_csv("lm_case3.csv",show_col_types = FALSE))
plot_component_heatmap <- function(data, component_indices, highlight_df) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(ggnewscale)
  library(patchwork)
  
  component_cols <- colnames(data)[component_indices]
  data[component_cols] <- lapply(data[component_cols], function(x) as.numeric(as.character(x)))
  
  # Merge highlight data
  data <- left_join(data, highlight_df, by = "metabolite")
  
  # Identify highlight status
  data <- data %>%
    mutate(
      Highlight = case_when(
        Beta_Dose != 0 & Beta_Serum != 0 ~ "Both",
        Beta_Dose != 0 ~ "Beta_Dose",
        Beta_Serum != 0 ~ "Beta_Serum",
        TRUE ~ "None"
      )
    )
  
  melted_data <- melt(data, id.vars = c("SUPER_PATHWAY", "SUB_PATHWAY", "metabolite", "Highlight"),
                      measure.vars = component_cols,
                      variable.name = "Component", value.name = "Loading")
  
  melted_data <- melted_data %>%
    mutate(
      RowID = factor(paste(SUPER_PATHWAY, SUB_PATHWAY, metabolite, sep = " | "),
                     levels = unique(paste(SUPER_PATHWAY, SUB_PATHWAY, metabolite, sep = " | "))),
      Component = factor(Component, levels = unique(Component)),
      x_pos = as.numeric(Component),
      x_super = -2,
      x_sub = -1,
      x_label = max(x_pos) + 1
    )
  
  super_colors <- scales::hue_pal()(length(unique(melted_data$SUPER_PATHWAY)))
  names(super_colors) <- unique(melted_data$SUPER_PATHWAY)
  
  sub_colors <- scales::hue_pal()(length(unique(melted_data$SUB_PATHWAY)))
  names(sub_colors) <- unique(melted_data$SUB_PATHWAY)
  
  text_colors <- c("None" = "black", "Beta_Dose" = "red", "Beta_Serum" = "green", "Both" = "purple")
  
  base_plot <- ggplot(melted_data, aes(y = RowID)) +
    geom_tile(aes(x = x_super, fill = SUPER_PATHWAY), width = 1) +
    scale_fill_manual(name = "SUPER_PATHWAY", values = super_colors) +
    ggnewscale::new_scale_fill() +
    
    geom_tile(aes(x = x_sub, fill = SUB_PATHWAY), width = 1) +
    scale_fill_manual(name = "SUB_PATHWAY", values = sub_colors) +
    ggnewscale::new_scale_fill() +
    
    geom_tile(aes(x = x_pos, fill = Loading), color = "white") +
    scale_fill_gradient2(name = "Loading", low = "blue", mid = "grey", high = "brown",
                         midpoint = 0, na.value = "grey",
                         limits = range(melted_data$Loading, na.rm = TRUE)) +
    
    geom_text(aes(x = x_label, label = metabolite, color = Highlight), hjust = 0, size = 3) +
    scale_color_manual(name = "Association with Drug", values = text_colors, na.translate = FALSE) +
    
    scale_x_continuous(
      breaks = c(unique(melted_data$x_pos), -1, -2),
      labels = c(levels(melted_data$Component), "SUB_PATHWAY", "SUPER_PATHWAY"),
      expand = expansion(mult = c(0.01, 0.6))
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(base_plot)
}
heatmap <- plot_component_heatmap(data=results,component_indices = 3:7,highlight_df = drug_data)
ggsave(plot = heatmap,filename = "heatmap.png",width=18,height=10,dpi = 600)
