library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
results <- read.csv("drug_non_drug_base_stats.csv")
results <- results %>% select(GLM_Beta,Logit_OR,DisplayName,Notes)
results$Logit_OR <- log10(results$Logit_OR)
results$GLM_Beta <- 10^(results$GLM_Beta)
results_metad <- read.csv("drug_non_drug_metadeconfound.csv")
plot_glm_vs_logit <- function(data,
                              x_col = "GLM_Beta",
                              y_col = "Logit_OR",
                              display_col = "DisplayName",
                              notes_col = "Notes") {
  # Extract second taxonomic level (phylum) from Notes
  data <- data %>%
    mutate(Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
           Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum))
  
  # Compute thresholds for top 1% of absolute values
  x_thresh <- quantile(abs(data[[x_col]]), 0.95, na.rm = TRUE)
  y_thresh <- quantile(abs(data[[y_col]]), 0.95, na.rm = TRUE)
  
  # Filter rows where either x or y is in the top 1%
  data_top1 <- data %>%
    filter(abs(!!sym(x_col)) >= x_thresh | abs(!!sym(y_col)) >= y_thresh)
  
  # Generate plot with color by Phylum
  ggplot(data_top1, aes_string(x = x_col, y = y_col, label = display_col, color = "Phylum")) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, color = "gray40", linetype = "dotted") +
    geom_text_repel(size = 3, max.overlaps = Inf) +
    labs(
      x = "10^GLM Beta Coefficient",
      y = "log10(Odds Ratio)",
      title = "Scatter Plot of Beta Coefficient from GLM and Logistic Regression",
      color = "Phylum"
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}
plot1 <- plot_glm_vs_logit(data=results)
plot_rho_bar <- function(data,
                         rho_col = "Rho",
                         display_col = "DisplayName",
                         notes_col = "Notes") {
  # Extract second taxonomic level (phylum) from Notes
  data <- data %>%
    mutate(Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
           Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum))
  
  # Compute threshold for top 5% of absolute Rho values
  rho_thresh <- quantile(abs(data[[rho_col]]), 0.7, na.rm = TRUE)
  
  # Filter top 5% entries
  data_top5 <- data %>%
    filter(abs(.data[[rho_col]]) >= rho_thresh)
  
  # Reorder DisplayName based on Rho for plotting
  data_top5 <- data_top5 %>%
    arrange(.data[[rho_col]]) %>%
    mutate(!!display_col := factor(.data[[display_col]], levels = .data[[display_col]]))
  
  # Create horizontal bar plot
  ggplot(data_top5, aes(x = .data[[rho_col]], y = .data[[display_col]], fill = Phylum)) +
    geom_col() +
    labs(
      x = "Rho",
      y = "CAGID",
      title = "Correlation Coefficient Plot (Top 30%)",
      fill = "Phylum"
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}
plot2 <- plot_rho_bar(data=results_metad)
