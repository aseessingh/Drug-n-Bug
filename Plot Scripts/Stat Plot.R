library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
results <- read.csv("drug_non_drug_base_mgs_stats.csv")
results <- results %>% select(GLM_Beta,Logit_OR,DisplayName,Notes,Rank)
results$Logit_OR <- log10(results$Logit_OR)
results$GLM_Beta <- 100*(results$GLM_Beta)
plot_glm_vs_logit <- function(data,
                              x_col = "GLM_Beta",
                              y_col = "Logit_OR",
                              label_col = "DisplayName",
                              fallback_label_col = "Species",
                              notes_col = "Notes",
                              group_col = "Phylum") {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggrepel)
  
  # Step 1: Filter for Rank == 1
  data <- data %>% filter(Rank == 1)
  
  # Step 2: Extract Phylum and Species from Notes
  data <- data %>%
    mutate(
      Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
      Species = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 7],
      Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum),
      Species = ifelse(is.na(Species) | Species == "", !!sym(label_col), Species),
      Label = ifelse(!is.na(!!sym(label_col)) & !!sym(label_col) != "", 
                     paste0(!!sym(label_col), "\n", Species), 
                     !!sym(fallback_label_col))
    )
  
  # Step 3: Compute thresholds for top 10%
  x_thresh <- quantile(abs(data[[x_col]]), 0.9, na.rm = TRUE)
  y_thresh <- quantile(abs(data[[y_col]]), 0.9, na.rm = TRUE)
  
  data_top <- data %>%
    filter(abs(!!sym(x_col)) >= x_thresh | abs(!!sym(y_col)) >= y_thresh)
  
  # Step 4: Define aesthetic mapping
  aesthetic_mapping <- aes_string(x = x_col, y = y_col)
  if (group_col == "Phylum") {
    aesthetic_mapping <- modifyList(aesthetic_mapping, aes_string(color = group_col))
  } else if (group_col == "Ontology") {
    aesthetic_mapping <- modifyList(aesthetic_mapping, aes_string(shape = group_col))
  }
  
  # Step 5: Plot
  ggplot(data_top, aesthetic_mapping) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, color = "gray40", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
    geom_text_repel(
      data = data_top %>% filter(!!sym(x_col) > 5 | !!sym(x_col) < -5),
      aes_string(label = "Label"),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      x = "GLM Beta Coefficient (as percentage, adjusted for demographics)",
      y = "log10(Odds Ratio) (adjusted for demographics)",
      title = "Beta Coefficient vs Odds Ratio Plot (Top 10%)",
      color = if (group_col == "Phylum") "Phylum" else NULL,
      shape = if (group_col == "Ontology") "Ontology" else NULL
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    expand_limits(x = 0, y = 0)  # Ensures (0, 0) is centered
}
plot1 <- plot_glm_vs_logit(data=results)

