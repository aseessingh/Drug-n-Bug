library(readr)
library(dplyr)
build_metabosignal_network <- function(compound_df, metabo_paths,signalling_paths, gene_name, organism_code = "hsa", n_cores = 8) {
  if (!requireNamespace("MetaboSignal", quietly = TRUE)) stop("Please install the MetaboSignal package.")
  if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("Please install the KEGGREST package.")
  
  library(MetaboSignal)
  library(parallel)
  library(KEGGREST)
  
  # Step 0: Preprocess compound_df
  compound_df <- compound_df[!is.na(compound_df[,1]) & trimws(compound_df[,1]) != "", , drop = FALSE]
  compound_df[,1] <- trimws(compound_df[,1])
  compound_df[,1] <- sapply(compound_df[,1], function(x) {
    ids <- unlist(strsplit(x, ","))
    ids <- trimws(ids)
    ids <- ifelse(grepl("^cpd:", ids), ids, paste0("cpd:", ids))  # keep cpd: prefix
    paste(ids, collapse = ",")
  })
  
  # Step 1: Parse and group compound IDs
  compound_groups <- mclapply(1:nrow(compound_df), function(i) {
    unlist(strsplit(compound_df[i, 1], ","))
  }, mc.cores = n_cores)
  names(compound_groups) <- rownames(compound_df)
  
  all_compounds <- unique(unlist(compound_groups))
  
  # Step 2: Build network using all pathways
  MS_table <- MS_keggNetwork(metabo_paths = metabo_paths, signaling_paths = signalling_paths)
  
  # Step 3: Merge isomers
  MS_table <- Reduce(function(tbl, group) {
    main_id <- group[1]
    if (length(group) > 1) {
      tbl <- MS_replaceNode(node1 = group[-1], node2 = main_id, tbl)
    }
    return(tbl)
  }, compound_groups, init = MS_table)
  
  # Step 4: Convert gene name to KEGG orthology ID
  gene_info <- MS_convertGene(genes = gene_name, organism_code = organism_code, organism_name = "human", output = "matrix")
  gene_kegg_id <- gene_info[1, "KEGG_ID"]
  
  # Step 5: Compute shortest paths
  target_compounds <- sapply(compound_groups, function(group) group[1])
  mapped_check <- MS_findMappedNodes(nodes = target_compounds, MS_table)
  if (length(mapped_check$mapped_nodes) == 0) {
    stop("❌ None of the target metabolites were mapped onto the network.")
  }
  
  distances <- MS_distances(MS_table, organism_code = organism_code,
                            source_genes = gene_kegg_id,
                            target_metabolites = target_compounds,
                            names = TRUE)
  
  # Step 6: Export network
  MS_shortestPathsNetwork(MS_table, organism_code = organism_code,
                          source_nodes = gene_kegg_id,
                          target_nodes = target_compounds,
                          type = "bw", file_name = "MetaboSignal_Network_2", names = TRUE)
  
  # Step 7: Create node label file with compound names
  original_names <- rownames(compound_df)
  main_ids <- sapply(compound_groups, function(group) group[1])
  name_to_kegg <- setNames(main_ids, original_names)
  
  # Fetch compound names from KEGG
  compound_names <- sapply(main_ids, function(cid) {
    tryCatch({
      entry <- keggGet(cid)[[1]]
      if (!is.null(entry$NAME)) {
        strsplit(entry$NAME[1], ";")[[1]][1]
      } else {
        gsub("^cpd:", "", cid)
      }
    }, error = function(e) {
      gsub("^cpd:", "", cid)
    })
  }, USE.NAMES = FALSE)
  
  node_labels <- data.frame(
    Node = gsub("^cpd:", "", main_ids),
    Label = compound_names,
    stringsAsFactors = FALSE
  )
  write.table(node_labels, file = "MetaboSignal_NodeLabels_2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(list(
    distances = distances,
    gene_kegg_id = gene_kegg_id,
    target_compounds = target_compounds,
    name_to_kegg = name_to_kegg,
    compound_names = compound_names
  ))
}
annotate_compounds_to_kegg <- function(df, column = "Compound", organism = "hsa", mc_cores = 2) {
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(KEGGREST)
  library(stringr)
  
  # Step 1: Clean and split compound IDs
  cleaned_compounds <- df %>%
    mutate(Compound = strsplit(as.character(.data[[column]]), ",")) %>%
    unnest(Compound) %>%
    mutate(Compound = str_trim(Compound)) %>%
    filter(Compound != "") %>%
    distinct()
  
  compound_ids <- cleaned_compounds$Compound
  
  if (length(compound_ids) == 0) {
    stop("No valid compound IDs found after cleaning.")
  }
  
  # Step 2: Get all organism-specific pathways
  org_pathways <- keggList("pathway", organism)
  org_pathways_df <- data.frame(
    PathwayID = names(org_pathways),
    PathwayName = unname(org_pathways),
    stringsAsFactors = FALSE
  )
  
  # Step 3: Define mapping function (no intersect)
  map_one_compound <- function(cid) {
    tryCatch({
      kegg_id <- paste0("cpd:", cid)
      generic_paths <- keggLink("pathway", kegg_id)
      generic_paths <- unname(generic_paths)
      
      data.frame(
        Compound = cid,
        PathwayID = generic_paths,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(Compound = cid, PathwayID = NA)
    })
  }
  
  # Step 4: Run in parallel
  results <- mclapply(compound_ids, map_one_compound, mc.cores = mc_cores)
  
  # Step 5: Combine results
  compound_pathways_df <- bind_rows(results)
  
  # Step 6: Return both data frames
  return(list(
    compound_pathways = compound_pathways_df,
    organism_pathways = org_pathways_df
  ))
}
mapping <- annotate_compounds_to_kegg(compound_ids,"Final.KEGG",mc_cores = 8,organism = "hsa")
mapping$compound_pathways$PathwayID <- gsub("path:map","hsa",mapping$compound_pathways$PathwayID)
mapping$linked_pathways <- intersect(mapping$compound_pathways$PathwayID,mapping$organism_pathways$PathwayID)
signalling_paths <- intersect(mapping$linked_pathways,c("hsa02010", "hsa02060", "hsa03070", "hsa02020", "hsa04010", "hsa04013", "hsa04016", "hsa04011", "hsa04012", 
                                                        "hsa04014", "hsa04015", "hsa04310", "hsa04330", "hsa04340", "hsa04341", "hsa04350", "hsa04390", "hsa04391", 
                                                        "hsa04392", "hsa04370", "hsa04371", "hsa04630", "hsa04064", "hsa04668", "hsa04066", "hsa04068", "hsa04020", 
                                                        "hsa04070", "hsa04072", "hsa04071", "hsa04024", "hsa04022", "hsa04151", "hsa04152", "hsa04150", "hsa04075", 
                                                        "hsa04080", "hsa04082", "hsa04081", "hsa04060", "hsa04061", "hsa04512", "hsa04514"))
metabo_paths <- intersect(mapping$linked_pathways,c("hsa00010", "hsa00020", "hsa00030", "hsa00040", "hsa00051", "hsa00052", "hsa00053", "hsa00061", "hsa00062", 
                                                    "hsa00071", "hsa00073", "hsa00074", "hsa00100", "hsa00120", "hsa00121", "hsa00130", "hsa00140", "hsa00190", 
                                                    "hsa00195", "hsa00196", "hsa00220", "hsa00230", "hsa00232", "hsa00240", "hsa00250", "hsa00253", "hsa00254", 
                                                    "hsa00260", "hsa00261", "hsa00270", "hsa00280", "hsa00290", "hsa00300", "hsa00310", "hsa00311", "hsa00320", 
                                                    "hsa00330", "hsa00331", "hsa00332", "hsa00333", "hsa00340", "hsa00350", "hsa00360", "hsa00363", "hsa00364", 
                                                    "hsa00365", "hsa00380", "hsa00400", "hsa00401", "hsa00402", "hsa00403", "hsa00404", "hsa00405", "hsa00410", 
                                                    "hsa00430", "hsa00440", "hsa00450", "hsa00460", "hsa00470", "hsa00480", "hsa00500", "hsa00510", "hsa00511", 
                                                    "hsa00512", "hsa00513", "hsa00514", "hsa00515", "hsa00520", "hsa00521", "hsa00522", "hsa00523", "hsa00524", 
                                                    "hsa00525", "hsa00531", "hsa00532", "hsa00533", "hsa00534", "hsa00540", "hsa00541", "hsa00542", "hsa00543", 
                                                    "hsa00550", "hsa00552", "hsa00561", "hsa00562", "hsa00563", "hsa00564", "hsa00565", "hsa00571", "hsa00572", 
                                                    "hsa00600", "hsa00601", "hsa00603", "hsa00604", "hsa00620", "hsa00621", "hsa00622", "hsa00623", "hsa00624", 
                                                    "hsa00625", "hsa00626", "hsa00627", "hsa00630", "hsa00633", "hsa00640", "hsa00642", "hsa00643", "hsa00650", 
                                                    "hsa00660", "hsa00670", "hsa00680", "hsa00710", "hsa00720", "hsa00730", "hsa00740", "hsa00750", "hsa00760", 
                                                    "hsa00770", "hsa00780", "hsa00785", "hsa00790", "hsa00791", "hsa00830", "hsa00860", "hsa00900", "hsa00901", 
                                                    "hsa00902", "hsa00903", "hsa00904", "hsa00905", "hsa00906", "hsa00907", "hsa00908", "hsa00909", "hsa00910", 
                                                    "hsa00920", "hsa00930", "hsa00940", "hsa00941", "hsa00942", "hsa00943", "hsa00944", "hsa00945", "hsa00946", 
                                                    "hsa00950", "hsa00960", "hsa00965", "hsa00966", "hsa00980", "hsa00981", "hsa00982", "hsa00983", "hsa00984", 
                                                    "hsa00996", "hsa00997", "hsa00998", "hsa00999", "hsa01010", "hsa01040", "hsa01051", "hsa01052", "hsa01053", 
                                                    "hsa01054", "hsa01055", "hsa01056", "hsa01057", "hsa01059", "hsa01060", "hsa01061", "hsa01062", "hsa01063", 
                                                    "hsa01064", "hsa01065", "hsa01066", "hsa01070", "hsa01100", "hsa01110", "hsa01120", "hsa01200", "hsa01210", 
                                                    "hsa01212", "hsa01220", "hsa01230", "hsa01232", "hsa01240", "hsa01250", "hsa01310", "hsa01320")
)
metabosignal <- build_metabosignal_network(compound_ids,metabo_paths,signalling_paths,"GLP1R")
metabosignal_GLP1R <- build_metabosignal_network(targets_GLP1R_try,metabo_paths,signalling_paths,"GLP1R")
compound_ids <- as.data.frame(read.csv("metabolites all kegg.csv"))
compound_ids <- na.omit(compound_ids)
targets_PRKAB1 <- as.data.frame(read_tsv("MetaboSignal_Network_TargetNodes.txt"))
targets_PRKAB1 <- targets_PRKAB1 %>% slice(-1) %>% select(node)
targets_PRKAB1$node <- gsub("cpd:","",targets_PRKAB1$node)
targets_GLP1R_try <- compound_ids %>%
  filter(!Final.KEGG %in% targets_PRKAB1$node)
replace_kegg_compounds <- function(input_file, output_file = "replaced_output.txt") {
  # Read the input file
  text <- readLines(input_file)
  
  # Extract all cpd:xyz patterns
  compound_ids <- str_extract_all(text, "cpd:[A-Za-z0-9]+") %>%
    unlist() %>%
    unique()
  
  # Query KEGG for each compound
  compound_names <- sapply(compound_ids, function(id) {
    tryCatch({
      entry <- keggGet(id)[[1]]
      entry$NAME[1]
    }, error = function(e) {
      NA
    })
  })
  
  # Replace each cpd:xyz with its name
  for (i in seq_along(compound_ids)) {
    if (!is.na(compound_names[i])) {
      text <- gsub(compound_ids[i], compound_names[i], text, fixed = TRUE)
    }
  }
  
  # Write the modified text to a new file
  writeLines(text, output_file)
  
  message("Replacement complete. Output saved to: ", output_file)
}
replace_kegg_compounds(input_file = "MetaboSignal_Network_Network.txt")
