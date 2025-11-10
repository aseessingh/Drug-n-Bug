# Load required packages
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}
library(KEGGREST)
library(readr)

# List of file paths
file_paths <- c(
  "drug_non_drug_base_KEGG_KO_stats.csv",
  "drug_non_drug_0_base_KEGG_KO_stats.csv",
  "drug_non_drug_base_KEGG_pathway_stats.csv",
  "drug_non_drug_0_base_KEGG_pathway_stats.csv",
  "drug_non_drug_base_KEGG_module_stats.csv",
  "drug_non_drug_0_base_KEGG_module_stats.csv"
)

# Initialize a list to store DisplayName values
all_display_names <- list()

# Read each file, filter Rank == 1, and extract DisplayName
for (file in file_paths) {
  if (file.exists(file)) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Filter for Rank == 1 if column exists
    if ("Rank" %in% colnames(df)) {
      df <- df[df$Rank == 1, ]
    }
    
    # Extract DisplayName if present
    if ("DisplayName" %in% colnames(df)) {
      all_display_names[[file]] <- df$DisplayName
    } else {
      warning(paste("DisplayName column not found in", file))
    }
  } else {
    warning(paste("File not found:", file))
  }
}

# Combine and find duplicates
combined_display_names <- unlist(all_display_names)
duplicate_names <- unique(combined_display_names[duplicated(combined_display_names)])
kegg <- data.frame(duplicate_names, stringsAsFactors = FALSE)

# Define row indices for each type (adjust as needed)
ko_rows <- c(1:439)
module_rows <- c(461:513)
pathway_rows <- c(440:460)

# Helper: Get KEGG name
get_kegg_name <- function(id) {
  tryCatch({
    entry <- keggGet(id)[[1]]
    if (!is.null(entry$NAME)) entry$NAME else "unknown"
  }, error = function(e) "not found")
}

# Helper: Get linked IDs
get_kegg_links <- function(target, source_id) {
  tryCatch({
    links <- keggLink(target, source_id)
    unique(sub(paste0(target, ":"), "", links))
  }, error = function(e) NA)
}

# Process each ID
results <- lapply(seq_len(nrow(kegg)), function(i) {
  id <- kegg$duplicate_names[i]
  if (i %in% ko_rows) {
    ko_id <- paste0("ko:", id)
    list(
      Input_ID = id,
      Type = "KO",
      Name = get_kegg_name(ko_id),
      Modules = paste(get_kegg_links("module", ko_id), collapse = "; "),
      Pathways = paste(get_kegg_links("pathway", ko_id), collapse = "; ")
    )
  } else if (i %in% module_rows) {
    module_id <- paste0("md:", id)
    list(
      Input_ID = id,
      Type = "Module",
      Name = get_kegg_name(module_id),
      Modules = id,
      Pathways = paste(get_kegg_links("pathway", module_id), collapse = "; ")
    )
  } else if (i %in% pathway_rows) {
    pathway_id <- paste0("path:", id)
    list(
      Input_ID = id,
      Type = "Pathway",
      Name = get_kegg_name(pathway_id),
      Modules = NA,
      Pathways = id
    )
  } else {
    list(
      Input_ID = id,
      Type = "Unknown",
      Name = "unknown",
      Modules = NA,
      Pathways = NA
    )
  }
})

# Convert to data frame and save
kegg_results <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
write.csv(kegg_results, "kegg_hierarchy_output.csv", row.names = FALSE)
