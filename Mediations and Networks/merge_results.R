library(dplyr)
library(parallel)

# Define the directory containing the result files
result_dir <- "Chunked Results"

# Define datasets and methods
datasets <- c("M1", "M3")
methods <- c("causal")

# Initialize results list
results <- setNames(
  lapply(datasets, function(ds) setNames(vector("list", length(methods)), methods)),
  datasets
)

# List all CSV files in the result directory
files <- list.files(result_dir, pattern = "\\.csv$", full.names = TRUE)

# Process each file
for (file in files) {
  for (dataset in datasets) {
    for (method in methods) {
      if (grepl(paste0("mediation_", dataset, "_chunk_\\d+_", method, "\\.csv$"), file)) {
        if (file.info(file)$size > 0) {
          df <- tryCatch(read.csv(file), error = function(e) NULL)
          if (!is.null(df)) {
            results[[dataset]][[method]] <- c(results[[dataset]][[method]], list(df))
          } else {
            message(paste("Skipping unreadable file:", file))
          }
        } else {
          message(paste("Skipping empty file:", file))
        }
      }
    }
  }
}

# Function to export combined and significant (acme_p_adj < 0.1) data
split_and_export_sig <- function(data, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  sig <- data %>% filter(acme_p_adj < 0.1)
  
  write.csv(data, file.path(output_dir, "combined.csv"), row.names = FALSE)
  write.csv(sig, file.path(output_dir, "sig.csv"), row.names = FALSE)
}

# Parallel processing using mclapply with return messages
mclapply(names(results), function(dataset) {
  lapply(names(results[[dataset]]), function(method) {
    if (length(results[[dataset]][[method]]) > 0) {
      combined <- bind_rows(results[[dataset]][[method]])
      output_dir <- file.path("results", dataset, method)
      split_and_export_sig(combined, output_dir)
      return(paste("Processed", dataset, method))
    } else {
      return(paste("No data for", dataset, method))
    }
  })
}, mc.cores = detectCores())

