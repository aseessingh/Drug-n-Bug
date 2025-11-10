library(dplyr)
indirect <- as.data.frame(read.csv("mediation group adjusted sig.csv"))
direct <- as.data.frame(read.csv("simple_output_significant_linear.csv"))
indirect <- indirect %>% filter(bug %in% bugs_of_interest$y_var)
direct <- direct %>% filter(microbe %in% bugs_of_interest$y_var)
variable_annotation <- as.data.frame(read.csv("hub.taxon_adjusted.MGS.down.10000000.variables.v6.r.csv"))
variable_annotation <- variable_annotation %>% select(VariableID,DisplayName)
metab_MGS_corr <- as.data.frame(read.csv("metab_MGS_corr.csv"))
build_and_plot_igraph_microbe_metabolite_network <- function(
    direct, indirect, correlation_matrix = NULL, display_key = NULL,
    export_prefix = "microbe_metabolite_network",
    layout_type = "fr",
    layout_scale = 1.5,
    jitter = FALSE,
    jitter_sd = 0.05,
    layout_charge = 0.05,
    layout_niter = 1000
) {
  library(dplyr)
  library(tidyr)
  library(igraph)
  
  # Prepare correlation matrix with symmetric lookup
  if (!is.null(correlation_matrix)) {
    correlation_matrix <- correlation_matrix %>%
      rename(microbe1 = 1, microbe2 = 2, cor_value = 5) %>%
      mutate(direction = ifelse(cor_value > 0, "Positive", "Negative"))
    
    correlation_matrix_rev <- correlation_matrix %>%
      transmute(microbe1 = microbe2, microbe2 = microbe1, cor_value, direction)
    
    correlation_matrix <- bind_rows(correlation_matrix, correlation_matrix_rev)
  }
  
  # Drug Related Microbe → Mediating Microbe
  edges_drug_to_microbe <- indirect %>%
    transmute(source = bug, target = microbe) %>%
    distinct()
  
  if (!is.null(correlation_matrix)) {
    edges_drug_to_microbe <- edges_drug_to_microbe %>%
      left_join(correlation_matrix, by = c("source" = "microbe1", "target" = "microbe2")) %>%
      mutate(direction = ifelse(is.na(direction), "Unknown", direction))
  } else {
    edges_drug_to_microbe <- edges_drug_to_microbe %>%
      mutate(direction = "Unknown")
  }
  
  edges_drug_to_microbe <- edges_drug_to_microbe %>%
    mutate(interaction = "Drug Related Microbe - Mediating Microbe")
  
  # Remove Unknown direction edges and downstream metabolite interactions
  valid_microbe_links <- edges_drug_to_microbe %>%
    filter(direction != "Unknown")
  
  valid_microbes <- unique(valid_microbe_links$target)
  
  edges_microbe_to_metabolite <- indirect %>%
    filter(microbe %in% valid_microbes) %>%
    transmute(source = microbe, target = metabolite,
              interaction = "Mediating Microbe - Metabolite",
              direction = ifelse(est_ACME > 0, "Positive", "Negative"))
  
  edges_drug_to_metabolite <- direct %>%
    transmute(source = microbe, target = metabolite,
              interaction = "Drug Related Microbe - Metabolite",
              direction = ifelse(beta > 0, "Positive", "Negative"))
  
  # Combine all edges
  edges <- bind_rows(valid_microbe_links, edges_microbe_to_metabolite, edges_drug_to_metabolite) %>%
    distinct() %>%
    filter(!grepl("^X-", source) & !grepl("^X-", target))
  
  # Create nodes
  nodes <- tibble(id = unique(c(edges$source, edges$target)))
  if (!is.null(display_key)) {
    nodes <- nodes %>%
      left_join(display_key, by = c("id" = "VariableID")) %>%
      mutate(label = ifelse(is.na(DisplayName), id, DisplayName))
  } else {
    nodes <- nodes %>% mutate(label = id)
  }
  
  nodes <- nodes %>%
    mutate(type = case_when(
      id %in% indirect$bug ~ "Drug Related Microbe",
      id %in% indirect$microbe ~ "Mediating Microbe",
      id %in% c(indirect$metabolite, direct$metabolite) ~ "Metabolite",
      TRUE ~ "Other"
    )) %>%
    filter(!grepl("^X-", id))
  
  # Build igraph object
  g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  
  # Set visual attributes
  V(g)$color <- case_when(
    V(g)$type == "Drug Related Microbe" ~ "yellow",
    V(g)$type == "Mediating Microbe" ~ "lightblue",
    V(g)$type == "Metabolite" ~ "orange",
    TRUE ~ "green"
  )
  
  V(g)$label <- V(g)$label
  deg <- degree(g, mode = "all")
  V(g)$size <- 5 + 3 * log(deg + 1)
  V(g)$label.cex <- 1.5
  V(g)$label.font <- 2
  
  E(g)$color <- case_when(
    E(g)$interaction == "Drug Related Microbe - Mediating Microbe" ~ "purple",
    E(g)$interaction == "Mediating Microbe - Metabolite" ~ "darkblue",
    E(g)$interaction == "Drug Related Microbe - Metabolite" ~ "pink",
    TRUE ~ "gray"
  )
  
  E(g)$lty <- ifelse(E(g)$direction == "Negative", 2, ifelse(E(g)$direction == "Positive", 1, 3))
  E(g)$arrow.size <- 0.6
  E(g)$width <- 2.5
  
  # Save plot
  png_file <- paste0(export_prefix, "_igraph_plot_30x18_600dpi.png")
  png(png_file, width = 30, height = 18, units = "in", res = 600)
  
  layout_func <- switch(layout_type,
                        fr = layout_with_fr,
                        circle = layout_in_circle,
                        kk = layout_with_kk,
                        lgl = layout_with_lgl,
                        graphopt = layout_with_graphopt,
                        drl = layout_with_drl,
                        mds = layout_with_mds,
                        gem = layout_with_gem,
                        random = layout_randomly,
                        star = layout_as_star,
                        layout_with_fr)
  layout_coords <- layout_func(g) * layout_scale
  
  if (jitter) {
    set.seed(123)
    layout_coords <- layout_coords + matrix(rnorm(length(layout_coords), 0, jitter_sd), ncol = 2)
  }
  
  plot(g,
       layout = layout_coords,
       vertex.label.color = "black",
       main = "Microbiome-Metabolite Network",
       vertex.label.cex = 1.5,
       cex.main = 3
  )
  
  legend("bottomleft",
         legend = c("Drug Related Microbe", "Mediating Microbe", "Metabolite", "Other"),
         col = c("yellow", "lightblue", "orange", "green"),
         pch = 21, pt.bg = c("yellow", "lightblue", "orange", "green"),
         pt.cex = 2.5, cex = 1.5, bty = "n", title = "Node Type")
  
  legend("bottomright",
         legend = c("Drug Related Microbe - Mediating Microbe", "Mediating Microbe - Metabolite", "Drug Related Microbe - Metabolite"),
         col = c("purple", "darkblue", "pink"),
         lty = 1, lwd = 3, cex = 1.5, bty = "n", title = "Edge Interaction")
  
  legend("topright",
         legend = c("Positive", "Negative", "Unknown"),
         lty = c(1, 2, 3), lwd = 3, col = "black", cex = 1.5, bty = "n", title = "Edge Direction")
  
  dev.off()
  
  # Create and export summary table (excluding correlation matrix columns)
  summary_table <- edges %>%
    select(source, target, interaction, direction) %>%
    left_join(nodes %>% select(id, source_type = type), by = c("source" = "id")) %>%
    left_join(nodes %>% select(id, target_type = type), by = c("target" = "id"))
  
  summary_file <- paste0(export_prefix, "_interaction_summary.csv")
  write.csv(summary_table, summary_file, row.names = FALSE)
  
  return(list(
    graph = g,
    plot_file = png_file,
    summary_file = summary_file
  ))
}
build_and_plot_igraph_microbe_metabolite_network(
    direct,
    indirect,
    display_key = variable_annotation,
    layout_type = "graphopt",layout_scale = 2000,jitter = TRUE,jitter_sd = 200,
    layout_charge = 0.2,layout_niter = 2000,correlation_matrix = metab_MGS_corr
  )

