# ---------------------------------------
#           Fuzzy Annotation Script
# ---------------------------------------

# --- Load required libraries ---
library(dplyr)
library(fuzzyjoin)
library(readr)

# --- INPUT: Set your file paths and parameters here ---
adj_file <- "hba1c_metabolites.csv"                       # Path to your 'adj' file (comma-separated)
chemical_annotation_file <- "chemical annotation.csv"          # Path to your annotation file (semicolon-separated)
max_distance <- 20                                              # Max string distance to attempt
method <- "osa"                                                 # String distance method: osa, jw, lv, etc.

# --- Read input data ---
adj <- read_csv(adj_file, locale = locale(encoding = "UTF-8")) %>%
  rename_with(trimws)

chemical_annotation <- read_csv(chemical_annotation_file, locale = locale(encoding = "UTF-8")) %>%
  rename_with(trimws)

# --- Function to perform fuzzy annotation ---
annotate_until_complete <- function(adj, chemical_annotation, max_distance = 20, method = "osa") {
  annotated_all <- tibble()
  remaining <- adj
  dist <- 0
  
  while (nrow(remaining) > 0 && dist <= max_distance) {
    message("Annotating with max_dist = ", dist)
    
    annotated_current <- stringdist_inner_join(
      remaining,
      chemical_annotation,
      by = c("...1" = "CHEMICAL_NAME"),
      method = method,
      max_dist = dist
    )
    
    if (nrow(annotated_current) == 0) {
      dist <- dist + 1
      next
    }
    
    annotated_all <- bind_rows(annotated_all, annotated_current) %>%
      distinct(`...1`, .keep_all = TRUE)
    
    remaining <- adj %>% filter(!(`...1` %in% annotated_all$`...1`))
    dist <- dist + 1
  }
  
  return(annotated_all)
}

# --- Run annotation ---
adj_annotated_final <- annotate_until_complete(adj, chemical_annotation, max_distance, method)
adj_annotated_final <- adj_annotated_final %>% select(CHEMICAL_NAME,SUPER_PATHWAY)
# --- Save result to CSV ---
write_csv(adj_annotated_final, "hba1c_metabolites_metadeconfoundR.csv")
