library(igraph)
library(dplyr)
library(readr)
extract_genus_from_notes <- function(notes_vector) {
  sapply(strsplit(notes_vector, ";"), function(x) {
    x <- trimws(x)
    x <- x[x != ""]
    
    if (length(x) >= 2) {
      genus <- x[length(x) - 1]
      
      # Remove "unclassified " prefix if present
      genus <- sub("^unclassified\\s+", "", genus, ignore.case = TRUE)
      
      # Check if genus is still empty or NA after cleanup
      if (is.na(genus) || genus == "" || grepl("unclassified", genus, ignore.case = TRUE)) {
        return("Unknown")
      } else {
        return(genus)
      }
    } else {
      return("Unknown")
    }
  })
}
microbe_annot <- read_csv("hub.taxon_adjusted.MGS.down.10000000.variables.v6.r.csv") %>%
  select(VariableID, Notes) %>%
  mutate(Genus = extract_genus_from_notes(Notes))
subpathway_df <- as.data.frame(read_csv("chemical annotation.csv"))
subpathway_df <- subpathway_df %>% select(CHEMICAL_NAME,SUB_PATHWAY)
build_correlation_network <- function(corr_file, subpathway_df = NULL, microbe_annot = NULL,
                                      combined_plot_file = "combined_network_plot.png",
                                      gephi_file = "correlation_network.graphml") {
  library(readr)
  library(dplyr)
  library(igraph)
  library(RColorBrewer)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(cowplot)
  library(ggforce)
  library(parallel)
  
  if (!file.exists(corr_file)) stop("Correlation file not found.")
  
  if (!is.null(subpathway_df) &&
      !all(c("CHEMICAL_NAME", "SUB_PATHWAY") %in% colnames(subpathway_df))) {
    stop("subpathway_df missing required columns.")
  }
  
  if (!is.null(microbe_annot) &&
      !all(c("VariableID", "Genus") %in% colnames(microbe_annot))) {
    stop("microbe_annot missing required columns.")
  }
  
  corr_df <- read_csv(corr_file, show_col_types = FALSE)
  
  # Precompute lookup vectors
  metabolite_map <- if (!is.null(subpathway_df)) setNames(subpathway_df$SUB_PATHWAY, subpathway_df$CHEMICAL_NAME) else character()
  microbe_map <- if (!is.null(microbe_annot)) setNames(microbe_annot$Genus, microbe_annot$VariableID) else character()
  
  annotate_vars <- function(df, var_col, type_col, group_col) {
    df %>%
      mutate(
        !!type_col := case_when(
          !!sym(var_col) %in% names(metabolite_map) ~ "metabolite",
          !!sym(var_col) %in% names(microbe_map) ~ "microbe",
          TRUE ~ NA_character_
        ),
        !!group_col := case_when(
          !!sym(var_col) %in% names(metabolite_map) ~ metabolite_map[!!sym(var_col)],
          !!sym(var_col) %in% names(microbe_map) ~ microbe_map[!!sym(var_col)],
          TRUE ~ NA_character_
        )
      )
  }
  
  # Parallel annotation
  n_cores <- detectCores(logical = FALSE)
  cl <- makeCluster(min(2, n_cores))
  clusterExport(cl, varlist = c("corr_df", "annotate_vars", "metabolite_map", "microbe_map"), envir = environment())
  clusterEvalQ(cl, library(dplyr))
  
  annotated_list <- parLapply(cl, list("var1", "var2"), function(var_col) {
    annotate_vars(corr_df, var_col, paste0(var_col, "_type"), paste0(var_col, "_group"))
  })
  stopCluster(cl)
  
  corr_df <- annotated_list[[1]]
  corr_df$var2_type <- annotated_list[[2]]$var2_type
  corr_df$var2_group <- annotated_list[[2]]$var2_group
  corr_df <- corr_df %>%
    mutate(edge_key = paste(pmin(var1, var2), pmax(var1, var2), sep = "_"))
  
  edges <- corr_df %>% select(var1, var2, corr, edge_key)
  nodes <- bind_rows(
    corr_df %>% select(name = var1, group = var1_group, type = var1_type),
    corr_df %>% select(name = var2, group = var2_group, type = var2_type)
  ) %>% distinct(name, .keep_all = TRUE)
  
  group_levels <- unique(na.omit(nodes$group))
  palette_size <- max(12, length(group_levels))
  group_colors <- setNames(colorRampPalette(brewer.pal(min(palette_size, 12), "Set3"))(length(group_levels)), group_levels)
  nodes$color <- rgb(t(col2rgb(group_colors[nodes$group])), maxColorValue = 255)
  
  g <- graph_from_data_frame(edges %>% select(var1, var2), vertices = nodes, directed = FALSE)
  V(g)$color <- nodes$color[match(V(g)$name, nodes$name)]
  V(g)$type <- nodes$type[match(V(g)$name, nodes$name)]
  V(g)$group <- nodes$group[match(V(g)$name, nodes$name)]
  
  edge_keys <- paste(pmin(as.character(get.edgelist(g)[, 1]), as.character(get.edgelist(g)[, 2])),
                     pmax(as.character(get.edgelist(g)[, 1]), as.character(get.edgelist(g)[, 2])), sep = "_")
  E(g)$weight <- abs(edges$corr[match(edge_keys, edges$edge_key)])
  E(g)$corr <- edges$corr[match(edge_keys, edges$edge_key)]
  E(g)$length <- 1 / abs(E(g)$corr + 1e-6)
  
  clusters <- cluster_louvain(g)
  V(g)$cluster <- as.character(membership(clusters)[V(g)$name])
  
  tg <- as_tbl_graph(g)
  layout <- create_layout(tg, layout = "fr", weights = E(g)$length)
  
  network_plot <- ggraph(layout) +
    geom_edge_link(aes(edge_width = weight), color = "grey60", alpha = 0.6) +
    geom_node_point(aes(color = group, shape = type), size = 3) +
    geom_mark_hull(aes(x = x, y = y, group = cluster, fill = cluster),
                   concavity = 5, alpha = 0.05, expand = unit(2, "mm")) +
    scale_shape_manual(values = c(microbe = 15, metabolite = 16), na.translate = FALSE) +
    scale_color_manual(values = group_colors, na.translate = FALSE) +
    theme_void() +
    theme(legend.position = "none")
  
  legend_df <- data.frame(group = names(group_colors), color = unname(group_colors)) %>%
    filter(!is.na(group))
  legend_plot <- ggplot(legend_df, aes(x = 1, y = reorder(group, desc(group)), color = group)) +
    geom_point(size = 5) +
    scale_color_manual(values = group_colors) +
    theme_void() +
    theme(legend.position = "none") +
    geom_text(aes(label = group), hjust = -0.1, color = "black") +
    xlim(1, 1.5)
  
  combined <- plot_grid(network_plot, legend_plot, rel_widths = c(4, 1), nrow = 1)
  ggsave(combined_plot_file, combined, width = 12, height = 8)
  
  write_graph(g, file = gephi_file, format = "graphml")
  
  return(list(graph = g, plot = combined))
}
network <- build_correlation_network(corr_file = c("metab_MGS_corr.csv"),subpathway_df = subpathway_df)
network_1 <- build_correlation_network(corr_file = "metab_KEGG_ko_corr.csv",combined_plot_file = "kegg_ko.png",gephi_file = "kegg_ko.graphml",subpathway_df = subpathway_df)
network_2 <- build_correlation_network(corr_file = "mgs_kegg_ko_corr.csv",combined_plot_file = "mgs_kegg_ko.png",gephi_file = "mgs_kegg_ko.graphml")
network_3 <- build_correlation_network(corr_file = c("metbugs_metagenome_corr.csv"),combined_plot_file = "metbugs_metagenome.png",gephi_file = "metbugs_metagenome.graphml")
kegg_int <- as.data.frame(read_csv("KEGG Entries clustering with metbugs.csv"))
kegg_annote <- as.data.frame(read_csv("hub.function_adjusted.KEGG_ko.down.10000000.fpkm.variables.v3.r.csv"))
kegg_hierarchy <- as.data.frame(read_csv("kegg_hierarchy_output.csv"))
kegg_annote <- kegg_annote %>% select(VariableID,DisplayName)
kegg_int <- kegg_int %>% left_join(kegg_annote,by=c("KEGG Entry"="VariableID")) %>% left_join(kegg_hierarchy,by=c("DisplayName"="Input_ID"))
write_csv(kegg_int,"KEGG Entries clustering with metbugs annote.csv")
