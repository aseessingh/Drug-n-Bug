library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
results <- read.csv("drug_non_drug_base_mgs_stats.csv")
results <- results %>% select(GLM_Beta,Logit_OR,Notes,DisplayName,Rank)
results$Logit_OR <- log10(results$Logit_OR)
results$GLM_Beta <- 100*(results$GLM_Beta)
plot_glm_vs_logit <- function(data,
                              x_col = "GLM_Beta",
                              y_col = "Logit_OR",
                              label_col = "Feature",
                              fallback_label_col = "Species",
                              notes_col = "Notes",
                              group_col = NULL,
                              drug_name = "Drug") {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggrepel)
  
  # Step 1: Filter for Rank == 1
  data <- data %>% filter(Rank == 1)
  
  # Step 2: Conditionally extract Phylum and Species from Notes
  if (notes_col %in% names(data)) {
    data <- data %>%
      mutate(
        Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
        Species = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 7],
        Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum),
        Species = ifelse(is.na(Species) | Species == "", !!sym(label_col), Species)
      )
    default_group_col <- "Phylum"
  } else {
    data <- data %>%
      mutate(
        Phylum = ifelse("Phylum" %in% names(data), Phylum, NA),
        Species = ifelse("Species" %in% names(data), Species, !!sym(label_col))
      )
    default_group_col <- if ("ontology" %in% names(data)) "ontology" else if ("Ontology" %in% names(data)) "Ontology" else NULL
  }
  
  # Step 3: Normalize ontology column to lowercase if it exists
  if ("ontology" %in% names(data)) {
    data <- data %>% mutate(ontology = tolower(ontology))
  } else if ("Ontology" %in% names(data)) {
    data <- data %>% mutate(Ontology = tolower(Ontology))
  }
  
  # Step 4: Create label (simplified to just Feature)
  data <- data %>%
    mutate(Label = !!sym(label_col))
  
  # Step 5: Compute thresholds for top 5%
  x_thresh <- quantile(abs(data[[x_col]]), 0.9, na.rm = TRUE)
  y_thresh <- quantile(abs(data[[y_col]]), 0.9, na.rm = TRUE)
  
  # Step 6: Define aesthetic mapping
  final_group_col <- if (!is.null(group_col)) group_col else default_group_col
  aesthetic_mapping <- aes_string(x = x_col, y = y_col)
  if (!is.null(final_group_col) && final_group_col %in% names(data)) {
    aesthetic_mapping <- modifyList(aesthetic_mapping, aes_string(color = final_group_col))
  }
  
  # Step 7: Plot
  ggplot(data, aesthetic_mapping) +
    geom_jitter(size = 3, alpha = 0.6, width = 0.2, height = 0.2) +
    geom_hline(yintercept = 0, color = "gray40", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
    geom_text_repel(
      data = data %>% filter(abs(!!sym(x_col)) >= x_thresh | abs(!!sym(y_col)) >= y_thresh),
      aes_string(label = "Label"),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2,
      max.iter = 3000,
      segment.size = 0.5,
      segment.curvature = 0.1,
      segment.ncp = 3
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("Positively Associated with", drug_name),
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "gray30") +
    annotate("text", x = -Inf, y = -Inf, label = paste("Negatively Associated with", drug_name),
             hjust = -0.1, vjust = -0.5, size = 4, fontface = "italic", color = "gray30") +
    labs(
      x = "GLM Beta Coefficient (as percentage, adjusted for demographics)",
      y = "log10(Odds Ratio) (adjusted for demographics)",
      title = "Beta Coefficient vs Odds Ratio Plot (Top 10%)",
      color = final_group_col
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.text = element_text(size = 12),       # Increased legend text size
      legend.title = element_text(size = 14)       # Increased legend title size
    ) +
    expand_limits(x = 0, y = 0)
}

plot2 <- plot_glm_vs_logit(data=results,drug_name = "Metformin",label_col = "DisplayName",notes_col = "Notes")
ggsave("gmm.png",plot = plot3,width=15,height=12,dpi=600)
