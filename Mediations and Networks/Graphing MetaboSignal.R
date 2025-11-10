plot_metabosignal_network_ggraph <- function(network_file, node_type_file, target_nodes_file, node_labels_file, layout_type = "sugiyama") {
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(readr)
  library(KEGGREST)
  
  supported_layouts <- c("fr", "kk", "lgl", "dh", "graphopt", "mds", "circle", "tree", "sugiyama", "grid", "randomly")
  if (!(layout_type %in% supported_layouts)) {
    stop(paste("Unsupported layout:", layout_type, "\nSupported layouts are:", paste(supported_layouts, collapse = ", ")))
  }
  
  network <- read_tsv(network_file, show_col_types = FALSE)
  node_types <- read_tsv(node_type_file, show_col_types = FALSE)
  target_nodes <- read_tsv(target_nodes_file, show_col_types = FALSE)
  node_labels <- read_tsv(node_labels_file, show_col_types = FALSE)
  
  node_labels$Node <- paste0("cpd:", node_labels$Node)
  compound_label_map <- setNames(node_labels$Label, node_labels$Node)
  
  compound_nodes <- unique(c(network$source, network$target))
  compound_nodes <- compound_nodes[grepl("^cpd:", compound_nodes)]
  missing_compounds <- setdiff(compound_nodes, names(compound_label_map))
  if (length(missing_compounds) > 0) {
    fetched_labels <- sapply(missing_compounds, function(id) {
      tryCatch({
        entry <- keggGet(id)[[1]]
        if (!is.null(entry$NAME)) entry$NAME[1] else id
      }, error = function(e) id)
    })
    compound_label_map <- c(compound_label_map, fetched_labels)
  }
  
  network$source_label <- ifelse(network$source %in% names(compound_label_map),
                                 compound_label_map[network$source], network$source)
  network$target_label <- ifelse(network$target %in% names(compound_label_map),
                                 compound_label_map[network$target], network$target)
  
  filtered_network <- network %>%
    select(source, target, type) %>%
    mutate(label = case_when(
      type == "k_activation" ~ "Activation",
      type == "k_compound:irreversible" ~ "Irreversible Conversion",
      type == "k_compound:reversible" ~ "Reversible Conversion",
      type == "k_expression" ~ "Expression Increase",
      type == "k_expression/indirect-effect" ~ "Indirect Expression Increase",
      type == "k_indirect-effect" ~ "Indirect Action",
      type == "k_indirect-effect/repression" ~ "Indirect Repression",
      TRUE ~ ""
    ))
  
  all_nodes <- unique(c(filtered_network$source, filtered_network$target))
  node_df <- tibble(name = all_nodes) %>%
    left_join(node_types, by = c("name" = "node")) %>%
    mutate(label = ifelse(name %in% names(compound_label_map),
                          compound_label_map[name], name),
           is_target = name %in% target_nodes$node)
  
  g <- tbl_graph(nodes = node_df, edges = filtered_network, directed = TRUE)
  
  p <- ggraph(g, layout = layout_type) +
    geom_edge_link(aes(color = label),
                   arrow = NULL,
                   end_cap = circle(3, 'mm')) +
    geom_node_point(aes(color = type), size = 4) +
    geom_node_text(aes(label = label, fontface = ifelse(is_target, "bold", "plain")),
                   repel = TRUE, size = 5, max.overlaps = Inf, box.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(name = "Type",
                       values = c("metabolic-gene" = "lightblue",
                                  "signaling-gene" = "lightgreen",
                                  "compound" = "black"),
                       labels = c("compound" = "Metabolite",
                                  "metabolic-gene" = "Metabolic Gene",
                                  "signaling-gene" = "Signalling Gene")) +
    theme_graph() +
    labs(title = "Host PRKAB1-Metabolite Interaction") +
    theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16)
    )
  
  output_file <- paste0("MetaboSignal_Network_", layout_type, ".png")
  ggsave(output_file, p, width = 30, height = 15, dpi = 600)
  
  return(output_file)
}
# Define all supported layouts
plot_metabosignal_network_ggraph(
    network_file = "MetaboSignal_Network_Network.txt",
    node_type_file = "MetaboSignal_Network_NodesType.txt",
    target_nodes_file = "MetaboSignal_Network_TargetNodes.txt",
    node_labels_file = "MetaboSignal_NodeLabels.txt",
    layout_type = "dh"
)

