library(dplyr)
library(readr)
library(parallel)
library(PMA)
library(VIM)

# Define file paths and corresponding variable names
file_paths <- list(
  metabolite = "MC_metabolites_wo_Drugs.csv",
  metadata = "metadata3.csv",
  mgs = "hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv",
  gmm = "hub.function_adjusted.GMM.down.10000000.fpkm.mgsamples.v5.r.csv",
  kegg_1 = "hub.function_adjusted.KEGG_ko.down.10000000.fpkm.mgsamples.v3.r.csv",
  kegg_2 = "hub.function_adjusted.KEGG_module.down.10000000.fpkm.mgsamples.v3.r.csv",
  kegg_3 = "hub.function_adjusted.KEGG_pathway.down.10000000.fpkm.mgsamples.v3.r.csv",
  group3 = "group 3.csv"
)

# Load all files in parallel using mclapply
ncores <- detectCores() - 1
data_list <- mclapply(file_paths, function(path) {
  if (grepl("MC_metabolites", path)) {
    as.data.frame(read_csv(path))
  } else {
    read.csv(path)
  }
}, mc.cores = ncores)

# Assign names and load into global environment
names(data_list) <- names(file_paths)
list2env(data_list, envir = .GlobalEnv)

# Preprocessing
group3 <- group3 %>% filter(SampleID %in% metadata$SampleID) %>% dplyr::select(SampleID, GLYCATHB)
metadata <- left_join(metadata, group3, by = "SampleID")
metabolite_AG <- metabolite %>% dplyr::select(MGS_PID, `1,5-anhydroglucitol (1,5-AG)`) %>% rename(AG = `1,5-anhydroglucitol (1,5-AG)`)
metadata <- left_join(metadata, metabolite_AG, by = c("SampleID" = "MGS_PID"))
metadata[, 7] <- kNN(metadata, variable = names(metadata)[7], k = 5)[, 7]
remove_M0_prefix <- function(df, column_name) {
  # Ensure column exists
  if (!column_name %in% names(df)) {
    stop("The specified column does not exist in the dataframe.")
  }
  
  # Remove "M0_" prefix only if it's at the start
  df[[column_name]] <- sub("^M0_", "", df[[column_name]])
  
  return(df)
}
gmm <- remove_M0_prefix(gmm,column_name = "MGSampleID")
kegg_1 <- remove_M0_prefix(kegg_1,column_name = "MGSampleID")
kegg_2 <- remove_M0_prefix(kegg_2,column_name = "MGSampleID")
kegg_3 <- remove_M0_prefix(kegg_3,column_name = "MGSampleID")
mgs <- remove_M0_prefix (mgs,column_name = "MGSampleID")

#Results File Loading
process_all_datasets <- function() {
  # Load required packages internally
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(parallel)
  })
  
  # Define the dataset types
  types <- c("GMM", "KEGG_KO", "KEGG_module", "KEGG_pathway", "metab", "mgs")
  
  # Define the processing function
  process_dataset <- function(type) {
    file1 <- paste0("drug_non_drug_0_base_", type, "_stats.csv")
    file2 <- paste0("drug_non_drug_base_", type, "_stats.csv")
    
    # Read and filter
    df1 <- tryCatch(read_csv(file1) %>% filter(Rank == 1), error = function(e) NULL)
    df2 <- tryCatch(read_csv(file2) %>% filter(Rank == 1), error = function(e) NULL)
    
    if (is.null(df1) || is.null(df2)) {
      warning(paste("Skipping", type, "due to read error"))
      return(NULL)
    }
    
    # Inner join by Feature
    joined <- inner_join(df1, df2, by = "Feature", suffix = c("_0_base", "_base"))
    
    # Add a column to indicate the source type
    joined <- joined %>% mutate(Source = type)
    
    return(joined)
  }
  
  # Use mclapply for parallel processing
  results_list <- mclapply(types, process_dataset, mc.cores = detectCores())
  names(results_list) <- types
  
  # Filter out NULLs and combine all into one data frame
  combined_results <- bind_rows(results_list[!sapply(results_list, is.null)])
  
  return(list(
    individual = results_list,
    combined = combined_results
  ))
}
# Run the function
joined_outputs <- process_all_datasets()
all_results_M1 <- as.data.frame(joined_outputs$combined)
all_data <- metabolite %>%
  dplyr::left_join(mgs,by=c("MGS_PID"="MGSampleID")) %>%
  dplyr::left_join(gmm,by=c("MGS_PID"="MGSampleID")) %>%
  dplyr::left_join(kegg_1,by=c("MGS_PID"="MGSampleID")) %>%
  dplyr::left_join(kegg_2,by=c("MGS_PID"="MGSampleID")) %>%
  dplyr::left_join(kegg_3,by=c("MGS_PID"="MGSampleID"))
all_data_M1 <- all_data %>% dplyr::select(MGS_PID,all_of(all_results_M1$Feature))
combined_df_M1 <- dplyr::left_join(metadata,all_data_M1,by=c("SampleID"="MGS_PID"))
convert_categorical_to_continuous <- function(df, col1_index, col2_index) {
  # Ensure col1_index is less than col2_index for correct reinsertion
  if (col1_index > col2_index) {
    temp <- col1_index
    col1_index <- col2_index
    col2_index <- temp
  }
  # Extract and one-hot encode the categorical columns
  cat1 <- model.matrix(~ . - 1, data = df[col1_index])
  cat2 <- model.matrix(~ . - 1, data = df[col2_index])
  
  # Remove original categorical columns
  df <- df[, -c(col1_index, col2_index)]
  
  # Reinsert the encoded columns at the correct positions
  df <- cbind(
    df[, 1:(col1_index - 1), drop = FALSE],
    cat1,
    df[, (col1_index):(col2_index - 2), drop = FALSE],
    cat2,
    df[, (col2_index - 1):ncol(df), drop = FALSE]
  )
  
  return(df)
}
combined_df_M1 <- convert_categorical_to_continuous(combined_df_M1,3,4)
combined_df_M1[,6] <- NULL
combined_df_M1[,9] <- NULL
combined_df_M1 <- combined_df_M1 %>% dplyr::select(SampleID,AGE,GENDERMale,GENDERFemale,CENTER_CDanemark,CENTER_CFrance,CENTER_CGermany,DOSAGE_METFORMIN_C,metformin,GLYCATHB,AG,everything())
combined_df_M1[,12:848] <- log1p(combined_df_M1[,12:848])
group_columns_by_category <- function(df) {
  # Define patterns for each category
  kegg_pattern <- "KEGG_ko"
  module_pattern <- "module"
  pathway_pattern <- "pathway"
  mgs_pattern <- "Taxon"
  
  # Identify columns by pattern
  kegg_cols <- grep(kegg_pattern, names(df), value = TRUE)
  module_cols <- grep(module_pattern, names(df), value = TRUE)
  pathway_cols <- grep(pathway_pattern, names(df), value = TRUE)
  taxon_cols <- grep(mgs_pattern, names(df), value = TRUE)
  
  # Assume remaining are metabolites
  all_cols <- names(df)
  metabolite_cols <- setdiff(all_cols, c(kegg_cols, module_cols, pathway_cols, taxon_cols))
  
  # Reorder columns: metabolites, KEGG, modules, pathways, taxon
  ordered_cols <- c(metabolite_cols, kegg_cols, module_cols, pathway_cols, taxon_cols)
  
  # Return reordered data frame
  df <- df[, ordered_cols]
  return(df)
}
combined_df_M1 <- group_columns_by_category(combined_df_M1)
write_csv(combined_df_M1, "combined_df_M1.csv")
#SparseCCA Dimensionality Reduction
select_mediators_scca_repeated_cv <- function(df, 
                                              x_indices, 
                                              y_index, 
                                              covar_indices, 
                                              mediator_indices, 
                                              penalty_grid, 
                                              nfolds = 10, 
                                              nrepeats = 3, 
                                              seed = 123) {
  set.seed(seed)
  
  # Extract matrices
  X <- scale(as.matrix(df[, x_indices]))
  Y <- scale(as.matrix(df[, y_index, drop = FALSE]))
  Covars <- scale(as.matrix(df[, covar_indices]))
  M <- scale(as.matrix(df[, mediator_indices]))
  
  # Residualize X and M with respect to covariates
  X_res <- residuals(lm(X ~ Covars))
  M_res <- residuals(lm(M ~ Covars))
  
  # Combine X and Y for canonical correlation
  XY <- cbind(X_res, Y)
  
  best_penalty <- NULL
  best_score <- -Inf
  
  for (penalty in penalty_grid) {
    scores <- numeric(nrepeats)
    
    for (r in 1:nrepeats) {
      folds <- sample(rep(1:nfolds, length.out = nrow(df)))
      fold_scores <- numeric(nfolds)
      
      for (f in 1:nfolds) {
        test_idx <- which(folds == f)
        train_idx <- setdiff(1:nrow(df), test_idx)
        
        cca_result <- tryCatch({
          CCA(x = XY[train_idx, ], z = M_res[train_idx, ], typex = "standard", typez = "standard", penaltyz = penalty)
        }, error = function(e) NULL)
        
        if (!is.null(cca_result)) {
          u_test <- XY[test_idx, ] %*% cca_result$u
          v_test <- M_res[test_idx, ] %*% cca_result$v
          fold_scores[f] <- cor(u_test, v_test)
        } else {
          fold_scores[f] <- NA
        }
      }
      
      scores[r] <- mean(fold_scores, na.rm = TRUE)
    }
    
    avg_score <- mean(scores, na.rm = TRUE)
    if (!is.na(avg_score) && avg_score > best_score) {
      best_score <- avg_score
      best_penalty <- penalty
    }
  }
  
  final_cca <- CCA(x = XY, z = M_res, typex = "standard", typez = "standard", penaltyz = best_penalty)
  selected_mediators <- which(abs(final_cca$v) > 1e-6)
  selected_names <- colnames(df)[mediator_indices[selected_mediators]]
  
  result_df <- df[, c(x_indices, y_index, covar_indices), drop = FALSE]
  if (length(selected_names) > 0) {
    result_df <- cbind(result_df, df[, selected_names, drop = FALSE])
  }
  attr(result_df, "penalty") <- best_penalty
  return(result_df)
}
select_mediators_scca_bootstrap <- function(df, 
                                            x_indices = 8:9, 
                                            y_index = 11, 
                                            covar_indices = 2:7, 
                                            mediator_indices,
                                            penalty_grid = c(0.01, 0.05, seq(0.1, 0.9, by = 0.01)), 
                                            nfolds = 10, 
                                            nrepeats = 3,
                                            nboot = 1000, 
                                            threshold = 0.3, 
                                            seed = 123,
                                            ncores = parallel::detectCores() - 1,
                                            output_csv = "selection_frequencies.csv",
                                            os_type = c("linux", "windows")) {
  library(PMA)
  library(dplyr)
  library(parallel)
  
  os_type <- match.arg(os_type)
  set.seed(seed)
  
  original_names <- colnames(df)[mediator_indices]
  cleaned_names <- make.names(original_names, unique = TRUE)
  name_key <- data.frame(original = original_names, cleaned = cleaned_names, stringsAsFactors = FALSE)
  colnames(df)[mediator_indices] <- cleaned_names
  
  mediator_names <- cleaned_names
  selection_matrix <- matrix(0, nrow = length(mediator_names), ncol = nboot)
  rownames(selection_matrix) <- name_key$original
  error_count <- 0
  penalty_used <- vector("numeric", nboot)
  
  bootstrap_function <- function(b) {
    set.seed(seed + b)
    boot_idx <- sample(1:nrow(df), replace = TRUE)
    df_boot <- df[boot_idx, ]
    
    tryCatch({
      result_df <- select_mediators_scca_repeated_cv(
        df = df_boot,
        x_indices = x_indices,
        y_index = y_index,
        covar_indices = covar_indices,
        mediator_indices = mediator_indices,
        penalty_grid = penalty_grid,
        nfolds = nfolds,
        nrepeats = nrepeats,
        seed = seed + b
      )
      
      selected_cleaned <- intersect(colnames(result_df), mediator_names)
      selected_original <- name_key$original[match(selected_cleaned, name_key$cleaned)]
      
      if ("penalty" %in% names(attributes(result_df))) {
        penalty_used[b] <<- attr(result_df, "penalty")
      }
      
      sel_vec <- rep(0, length(mediator_names))
      names(sel_vec) <- name_key$original
      sel_vec[selected_original] <- 1
      sel_vec
    }, error = function(e) {
      error_count <<- error_count + 1
      penalty_used[b] <<- NA
      rep(0, length(mediator_names))
    })
  }
  
  if (os_type == "linux") {
    results <- mclapply(1:nboot, mc.cores = ncores, FUN = bootstrap_function)
  } else {
    cl <- makeCluster(ncores)
    clusterExport(cl, varlist = c("df", "x_indices", "y_index", "covar_indices", "mediator_indices",
                                  "penalty_grid", "nfolds", "nrepeats", "seed", "mediator_names",
                                  "select_mediators_scca_repeated_cv", "name_key", "error_count", "penalty_used"), envir = environment())
    results <- parLapply(cl, 1:nboot, bootstrap_function)
    stopCluster(cl)
  }
  
  selection_matrix <- do.call(cbind, results)
  selection_freq <- rowMeans(selection_matrix)
  
  if (length(selection_freq) == 0 || is.null(names(selection_freq))) {
    warning("Selection frequency vector is empty or unnamed. No CSV will be written.")
    write.csv(data.frame(Mediator = character(), Selection_Frequency = numeric()), 
              file = output_csv, row.names = FALSE)
    return(df[, 1:12])
  }
  
  freq_df <- data.frame(Mediator = names(selection_freq), Selection_Frequency = selection_freq)
  write.csv(freq_df, file = output_csv, row.names = FALSE)
  
  stable_mediators <- names(selection_freq[selection_freq >= threshold])
  
  if (length(stable_mediators) < 0.1 * length(mediator_names)) {
    message("Fewer than 10% of mediators met the threshold. Applying fallback: selecting top 25% most frequent mediators.")
    n_fallback <- max(1, floor(0.25 * length(mediator_names)))
    top_mediators <- names(sort(selection_freq, decreasing = TRUE))[1:n_fallback]
    stable_mediators <- top_mediators
  }
  
  output_df <- df[, 1:12]
  if (length(stable_mediators) > 0) {
    cleaned_stable <- name_key$cleaned[match(stable_mediators, name_key$original)]
    output_df <- cbind(output_df, df[, cleaned_stable, drop = FALSE])
    colnames(output_df)[(ncol(output_df) - length(cleaned_stable) + 1):ncol(output_df)] <- stable_mediators
  }
  
  penalty_summary <- na.omit(penalty_used)
  if (length(penalty_summary) > 0) {
    most_common_penalty <- names(sort(table(penalty_summary), decreasing = TRUE))[1]
    message(paste("Most commonly selected penalty across bootstraps:", most_common_penalty))
  } else {
    message("No penalty information was recorded.")
  }
  
  message(paste("Bootstrap errors:", error_count, "out of", nboot))
  message(paste("Selection frequencies saved to:", output_csv))
  message(paste(length(stable_mediators), "mediators selected."))
  
  return(output_df)
}

combined_df_M3 <- select_mediators_scca_bootstrap(combined_df_M1,mediator_indices = 12:848,os_type="linux")
combined_df_M3 <- group_columns_by_category(combined_df_M3)
#combined_df_M3A <- select_mediators_scca(combined_df_M2A,mediator_indices = 12:1289)
#combined_df_M3B <- select_mediators_scca(combined_df_M2B,mediator_indices = 12:254)
write_csv(combined_df_M3,"combined_df_M3.csv")
#write_csv(combined_df_M3A,"combined_df_M3A.csv")
#write_csv(combined_df_M3B,"combined_df_M3B.csv")

