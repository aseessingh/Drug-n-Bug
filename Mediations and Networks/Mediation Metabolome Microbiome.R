library(readr)
library(dplyr)
library(parallel)
library(stringr)
library(tidyr)
metabolite_df <- as.data.frame(read_csv("metabolites.csv"))
metabolite_df <- metabolite_df %>% select(-c(1,3:8))
microbe_df <-  as.data.frame(read_csv("MGS.csv"))
microbe_df <- microbe_df %>% select(-1)
covar_df <- as.data.frame(read_csv("covar df.csv"))
covar_df <- covar_df %>% select(-1)
corr_df <- as.data.frame(read_csv("corr met MGS.csv",show_col_types = FALSE))
bugs_of_interest <- corr_df %>% filter(p_adj<0.05) %>% select(y_var)
bug_df <- microbe_df %>% select(1,bugs_of_interest$y_var)
linked_df <- read_csv("linked_microbes_combined.csv", show_col_types = FALSE) %>%
  mutate(connected_to = str_split(connected_to, "; ")) %>%
  unnest(connected_to) %>%
  rename(bug = connected_to)
run_mediation_analysis_lavaan <- function(
    metabolite_df,
    microbe_df,
    covar_df,
    linked_df,
    output_prefix = "mediation_results",
    n_boot = 1000,
    n_cores = parallel::detectCores()
) {
  suppressPackageStartupMessages({
    library(lavaan)
    library(dplyr)
    library(readr)
    library(purrr)
    library(tibble)
    library(stringr)
    library(parallel)
    library(tidyr)
    library(glue)
  })
  
  # Store original names and sanitize
  name_key <- tibble(
    original = unique(c(
      colnames(metabolite_df)[-1],
      colnames(microbe_df)[-1],
      colnames(covar_df)[-1],
      unlist(linked_df)
    ))
  ) %>%
    mutate(safe = make.names(original))
  
  sanitize_cols <- function(df) {
    colnames(df)[-1] <- name_key$safe[match(colnames(df)[-1], name_key$original)]
    df
  }
  
  metabolite_df <- sanitize_cols(metabolite_df %>% rename(sample_id = 1))
  microbe_df <- sanitize_cols(microbe_df %>% rename(sample_id = 1))
  covar_df <- sanitize_cols(covar_df %>% rename(sample_id = 1))
  linked_df <- linked_df %>% mutate(across(everything(), make.names))
  
  # Drop collinear covariates (correlation > 0.8)
  if (ncol(covar_df) > 2) {
    numeric_covars <- covar_df[, -1] %>% select(where(is.numeric))
    cor_matrix <- cor(numeric_covars, use = "pairwise.complete.obs")
    high_corr <- which(abs(cor_matrix) > 0.8 & lower.tri(cor_matrix), arr.ind = TRUE)
    to_drop <- unique(colnames(cor_matrix)[high_corr[, 2]])
    if (length(to_drop) > 0) {
      covar_df <- covar_df %>% select(-all_of(to_drop))
      message(glue("Dropped collinear covariates: {paste(to_drop, collapse = ', ')}"))
    }
  }
  
  full_data <- reduce(list(metabolite_df, microbe_df, covar_df), full_join, by = "sample_id")
  
  debug_log <- list()
  
  run_triplet <- function(i) {
    row <- linked_df[i, ]
    bug <- row$bug
    microbe <- row$microbe
    metabolite <- row$metabolite
    
    if (!all(c(bug, microbe, metabolite) %in% colnames(full_data))) {
      debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = "Missing variable(s)")
      return(NULL)
    }
    
    vars <- full_data[, c(bug, microbe, metabolite)]
    if (any(sapply(vars, function(x) all(is.na(x) | x == x[1])))) {
      debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = "Constant or missing variable(s)")
      return(NULL)
    }
    
    if (bug == microbe) {
      debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = "Bug and Microbe are the same")
      return(NULL)
    }
    
    covars <- colnames(covar_df)[-1]
    covar_str <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
    
    model <- paste0(
      microbe, " ~ a*", bug, if (covar_str != "1") paste0(" + ", covar_str) else "", "\n",
      metabolite, " ~ b*", microbe, " + cp*", bug, if (covar_str != "1") paste0(" + ", covar_str) else "", "\n",
      "ACME := a*b\n",
      "ADE := cp\n",
      "TotalEffect := cp + (a*b)\n",
      "PropMediated := ACME / TotalEffect"
    )
    
    fit <- tryCatch({
      sem(model, data = full_data, se = "bootstrap", bootstrap = n_boot)
    }, error = function(e) {
      debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = paste("SEM failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      result <- tryCatch({
        parameterEstimates(fit, standardized = TRUE) %>%
          filter(label %in% c("ACME", "ADE", "TotalEffect", "PropMediated")) %>%
          select(label, est, se, pvalue) %>%
          pivot_wider(names_from = label, values_from = c(est, se, pvalue)) %>%
          mutate(bug = bug, microbe = microbe, metabolite = metabolite)
      }, error = function(e) {
        debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = paste("Estimate extraction failed:", e$message))
        return(NULL)
      })
      return(result)
    } else {
      debug_log[[length(debug_log) + 1]] <<- tibble(i, bug, microbe, metabolite, reason = "Model did not converge")
      return(NULL)
    }
  }
  
  results_list <- mclapply(1:nrow(linked_df), run_triplet, mc.cores = n_cores)
  results_list <- results_list[!sapply(results_list, is.null)]
  results <- bind_rows(results_list)
  
  # Restore original names
  restore_names <- function(x) {
    name_key_map <- setNames(name_key$original, name_key$safe)
    cols_to_restore <- intersect(c("bug", "microbe", "metabolite"), colnames(x))
    if (length(cols_to_restore) > 0) {
      x <- x %>%
        mutate(across(all_of(cols_to_restore), ~ ifelse(. %in% names(name_key_map), name_key_map[.], .)))
    }
    return(x)
  }
  
  results <- restore_names(results)
  debug_df <- restore_names(bind_rows(debug_log))
  
  if ("pvalue_ACME" %in% colnames(results)) {
    results <- results %>%
      mutate(pvalue_ACME_fdr = p.adjust(pvalue_ACME, method = "fdr"))
  } else {
    warning("No valid mediation results with pvalue_ACME found.")
    results$pvalue_ACME_fdr <- NA
  }
  
  write_csv(results, paste0(output_prefix, "_all.csv"))
  sig_results <- results %>% filter(!is.na(pvalue_ACME_fdr) & pvalue_ACME_fdr < 0.05)
  write_csv(sig_results, paste0(output_prefix, "_significant.csv"))
  write_csv(debug_df, paste0(output_prefix, "_debug_log.csv"))
  
  return(list(all = results, significant = sig_results, debug_log = debug_df))
}
run_mediation_analysis_lavaan(metabolite_df,microbe_df,covar_df,linked_df,n_boot = 1000)


