library(MOFA2)
library(tidyverse)
library(readr)
library(lavaan)
library(parallel)
library(ppcor)

# Define datasets and omics ranges
datasets <- list(
  list(name = "M1", file = "combined_df_M1.csv", omics_ranges = list(
    GMM = 12:13,
    Metabolomics = 14:176,
    KEGG_KO = 177:615,
    KEGG_module = 616:668,
    KEGG_pathway = 669:689,
    MGS = 690:848
  )),
  list(name = "M3", file = "combined_df_M3.csv", omics_ranges = list(
    GMM = 12:12,
    Metabolomics = 13:128,
    KEGG_KO = 129:163,
    KEGG_module = 164:165,
    KEGG_pathway = 166:166,
    MGS = 167:221
  ))
)

# SEM column definitions
x_cols <- 8:9
y_cols <- 10:11
covar_cols <- 2:7

prepare_datasets <- function(datasets, x_cols, y_cols, covar_cols, threshold = 0.85) {
  for (i in seq_along(datasets)) {
    df <- readr::read_csv(datasets[[i]]$file)
    
    # Extract covariate columns
    covar_data <- df[, covar_cols]
    
    # Compute correlation matrix and identify highly correlated variables
    cor_matrix <- cor(covar_data, use = "pairwise.complete.obs")
    cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
    high_cor <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
    
    # Determine variables to drop
    vars_to_drop <- unique(colnames(cor_matrix)[high_cor[, 2]])
    
    # Save dropped covariates to CSV
    if (length(vars_to_drop) > 0) {
      drop_log <- data.frame(Dropped_Covariates = vars_to_drop)
      write.csv(drop_log, paste0("dropped_covariates_", datasets[[i]]$name, ".csv"), row.names = FALSE)
    }
    
    covar_names <- names(df)[covar_cols]
    cleaned_names <- setdiff(covar_names, vars_to_drop)
    cleaned_indices <- which(names(df) %in% cleaned_names)
    
    datasets[[i]]$data <- df
    datasets[[i]]$x_vars <- x_cols
    datasets[[i]]$y_vars <- y_cols
    datasets[[i]]$covar_vars <- cleaned_indices
  }
  return(datasets)
}
run_mofa_and_mediation <- function(ds) {
  library(MOFA2)
  library(lavaan)
  library(tidyverse)
  
  # Prepare MOFA object
  omics_ranges <- ds$omics_ranges
  df <- ds$data
  views <- lapply(omics_ranges, function(rng) t(as.matrix(df[, rng])))
  names(views) <- names(omics_ranges)
  
  mofa_obj <- create_mofa(views)
  model_opts <- get_default_model_options(mofa_obj)
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$maxiter <- 1000
  data_opts <- get_default_data_options(mofa_obj)
  
  mofa_obj <- prepare_mofa(
    mofa_obj,
    model_options = model_opts,
    training_options = train_opts,
    data_options = data_opts
  )
  
  mofa_obj <- run_mofa(mofa_obj, outfile = paste0("MOFA_", ds$name, ".hdf5"), use_basilisk = TRUE)
  
  # Save latent factors in long format
  factors <- get_factors(mofa_obj, factors = "all", as.data.frame = TRUE)
  write.csv(factors, paste0("latent_factors_", ds$name, ".csv"), row.names = FALSE)
  
  # Also save wide format for downstream use
  factors_wide <- factors %>%
    pivot_wider(names_from = factor, values_from = value)
  write.csv(factors_wide, paste0("latent_factors_", ds$name, "_wide.csv"), row.names = FALSE)
  
  # Save weights and variance explained
  weights <- get_weights(mofa_obj, views = "all", factors = "all", as.data.frame = TRUE)
  write.csv(weights, paste0("weights_", ds$name, ".csv"), row.names = FALSE)
  
  var_exp <- get_variance_explained(mofa_obj)
  var_df <- as.data.frame(var_exp$r2_per_factor)
  write.csv(var_df, paste0("variance_explained_", ds$name, ".csv"), row.names = TRUE)
  
  # Mediation analysis
  results <- list()
  factors_wide <- read.csv(paste0("latent_factors_", ds$name, "_wide.csv"))
  rownames(factors_wide) <- factors_wide$sample
  
  for (x in ds$x_vars) {
    for (y in ds$y_vars) {
      for (f in setdiff(colnames(factors_wide), "sample")) {
        sem_df <- cbind(df[, c(x, y, ds$covar_vars)], M = factors_wide[[f]])
        colnames(sem_df)[1:2] <- c("X", "Y")
        model <- '
          M ~ a*X
          Y ~ b*M + c*X
          indirect := a*b
          direct := c
          total := c + (a*b)
          proportion := (a*b)/(c + a*b)
        '
        fit <- tryCatch(sem(model, data = sem_df, se = "bootstrap", bootstrap = 1000), error = function(e) NULL)
        if (!is.null(fit)) {
          summary_df <- parameterEstimates(fit) %>%
            filter(label %in% c("indirect", "direct", "total", "proportion")) %>%
            dplyr::select(label, est, pvalue) %>%
            mutate(Factor = f, X = names(df)[x], Y = names(df)[y])
          results[[length(results) + 1]] <- summary_df
        }
      }
    }
  }
  
  final_df <- bind_rows(results) %>%
    group_by(X, Y) %>%
    mutate(p_adj = p.adjust(pvalue, method = "BH")) %>%
    ungroup()
  
  write.csv(final_df, paste0("mediation_summary_", ds$name, ".csv"), row.names = FALSE)
  sig_df <- final_df %>% filter(p_adj < 0.05)
  write.csv(sig_df, paste0("significant_mediations_", ds$name, ".csv"), row.names = FALSE)
}
run_partial_spearman <- function(ds) {
  library(ppcor)
  library(tidyverse)
  library(parallel)
  
  sem_data <- ds$data
  latent_factors <- read.csv(paste0("latent_factors_", ds$name, "_wide.csv"))
  rownames(latent_factors) <- latent_factors$sample
  latent_vars <- setdiff(colnames(latent_factors), "sample")
  y_vars <- ds$y_vars
  covar_vars <- ds$covar_vars
  
  combos <- expand.grid(Outcome = y_vars, Factor = latent_vars, stringsAsFactors = FALSE)
  
  cor_fun <- function(i) {
    yv <- combos$Outcome[i]
    f <- combos$Factor[i]
    test <- tryCatch({
      pcor.test(latent_factors[[f]], sem_data[[yv]], sem_data[, covar_vars], method = "spearman")
    }, error = function(e) list(estimate = NA, p.value = NA))
    data.frame(Factor = f, Outcome = yv, rho = test$estimate, p = test$p.value)
  }
  
  n_cores <- max(1, detectCores() - 1)
  cor_results_list <- mclapply(seq_len(nrow(combos)), cor_fun, mc.cores = n_cores)
  cor_df <- do.call(rbind, cor_results_list)
  
  cor_df <- cor_df %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    mutate(sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    mutate(label = paste0(sprintf("%.2f", as.numeric(rho)), sig))
  
  write.csv(cor_df, paste0("partial_spearman_correlations_", ds$name, ".csv"), row.names = FALSE)
}

# Step 0: Prepare datasets
datasets <- prepare_datasets(datasets, x_cols, y_cols, covar_cols)

# Step 1–3: Define wrapper
wrapper <- function(ds) {
  run_mofa_and_mediation(ds)
  run_partial_spearman(ds)
}

# Step 4: Run in parallel
n_cores <- detectCores()
mclapply(datasets, wrapper, mc.cores = n_cores)


