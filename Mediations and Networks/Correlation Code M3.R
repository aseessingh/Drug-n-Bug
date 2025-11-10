library(ppcor)
library(readr)
library(tidyverse)
library(parallel)
metagenome <- as.data.frame(read_csv("metagenome combined.csv"))
metagenome[,1] <- NULL
dataset <- as.data.frame(read_csv("combined_df_M1.csv"))
F7_MOFA <- as.data.frame(read_csv("M3_F7_TopQuart.csv"))
causalmed <- as.data.frame(read_csv("sig.csv"))
F7_MOFA <- F7_MOFA %>% select(feature)
causalmed <- causalmed %>% select(mediator) %>% distinct(mediator,.keep_all = TRUE)
causalmed <- causalmed[-c(1,47:53),]
causalmed <- as.data.frame(causalmed)
causalmed <- causalmed %>% rename(feature = causalmed)
met_mediators <- rbind(causalmed,F7_MOFA)
met_mediators <- met_mediators %>% distinct(feature,.keep_all = TRUE)
write.csv(met_mediators,"met mediators.csv")
dataset_met <- dataset %>% select(1:7,all_of(met_mediators$feature))
split_R_by_feature_type <- function(R) {
  substrings <- list(
    MGS = "Taxon",
    KEGG_ko = "KEGG_ko",
    KEGG_pathway = "KEGG_pathway",
    KEGG_module = "KEGG_module",
    GMM = "GMM"
  )
  
  split_dfs <- list()
  
  for (label in names(substrings)) {
    patterns <- substrings[[label]]
    
    # Combine multiple patterns into a single regex if needed
    if (length(patterns) > 1) {
      pattern <- paste(patterns, collapse = "|")
    } else {
      pattern <- patterns
    }
    
    # Find matching columns
    matching_cols <- grep(pattern, names(R), value = TRUE)
    
    if (length(matching_cols) > 0) {
      split_df <- R[, c("MGSampleID", matching_cols), drop = FALSE]
      split_dfs[[label]] <- split_df
    }
  }
  
  return(split_dfs)
}
combined <- split_R_by_feature_type(metagenome)
partial_spearman_corr <- function(x_df, y_df, covar_indices, x_indices, y_indices, n_cores = 8, thresh = 0.4, include_within = TRUE) {
  library(ppcor)
  library(parallel)
  library(dplyr)
  
  # Match by first column (assumed to be SampleID)
  sample_col_x <- x_df[[1]]
  sample_col_y <- y_df[[1]]
  common_samples <- intersect(sample_col_x, sample_col_y)
  
  x_df <- x_df[sample_col_x %in% common_samples, , drop = FALSE]
  y_df <- y_df[sample_col_y %in% common_samples, , drop = FALSE]
  
  # Reorder to match
  x_df <- x_df[match(common_samples, x_df[[1]]), ]
  y_df <- y_df[match(common_samples, y_df[[1]]), ]
  
  # Extract covariates and correlation variables
  covar_data <- x_df[, covar_indices, drop = FALSE]
  x_corr_df <- x_df[, x_indices, drop = FALSE]
  y_corr_df <- y_df[, y_indices, drop = FALSE]
  
  # Define group combinations
  if (include_within) {
    combos <- list(
      c("x", "x"),
      c("y", "y"),
      c("x", "y")  # only one direction
    )
  } else {
    combos <- list(
      c("x", "y")
    )
  }
  
  # Generate index combinations
  get_combos <- function(g1, g2) {
    d1 <- if (g1 == "x") x_corr_df else y_corr_df
    d2 <- if (g2 == "x") x_corr_df else y_corr_df
    expand.grid(i = seq_along(d1), j = seq_along(d2)) %>%
      filter(i < j | g1 != g2) %>%
      mutate(group1 = g1, group2 = g2)
  }
  
  all_combos <- do.call(rbind, Map(function(g) get_combos(g[1], g[2]), combos))
  
  # Correlation computation function
  compute_corr <- function(row) {
    g1 <- row$group1
    g2 <- row$group2
    i <- row$i
    j <- row$j
    
    d1 <- if (g1 == "x") x_corr_df else y_corr_df
    d2 <- if (g2 == "x") x_corr_df else y_corr_df
    
    x_vec <- d1[[i]]
    y_vec <- d2[[j]]
    
    data_mat <- data.frame(x = x_vec, y = y_vec, covar_data)
    
    result <- tryCatch({
      pcor.test(data_mat$x, data_mat$y, data_mat[, -c(1, 2)], method = "spearman")
    }, error = function(e) {
      return(list(estimate = NA, p.value = NA))
    })
    
    return(c(
      var1 = colnames(d1)[i],
      var2 = colnames(d2)[j],
      group1 = g1,
      group2 = g2,
      corr = result$estimate,
      p = result$p.value
    ))
  }
  
  # Parallel execution
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("x_corr_df", "y_corr_df", "covar_data", "compute_corr"), envir = environment())
    results <- parLapply(cl, split(all_combos, seq(nrow(all_combos))), compute_corr)
    stopCluster(cl)
  } else {
    results <- mclapply(split(all_combos, seq(nrow(all_combos))), compute_corr, mc.cores = n_cores)
  }
  
  # Format results
  results_df <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
  results_df$corr <- as.numeric(results_df$corr)
  results_df$p <- as.numeric(results_df$p)
  results_df$p_adj <- p.adjust(results_df$p, method = "fdr")
  
  filtered <- results_df %>%
    filter(p_adj < 0.05, abs(corr) > thresh)
  
  return(filtered)
}
metab_MGS_corr <- partial_spearman_corr(x_df = dataset,y_df = combined$MGS,x_indices = 8:91,y_indices = 2:210,covar_indices = 2:7,n_cores=8,thresh = 0.4,include_within = TRUE)
write_csv(metab_MGS_corr, "metab_MGS_corr.csv")
metagenome_wo_mgs <- cbind(combined$KEGG_ko,combined$KEGG_pathway,combined$KEGG_module,combined$GMM)
write.csv(metagenome_wo_mgs,"metagenome without cags.csv")

#Factor-Metabolites-MGS
factor_mofa <- as.data.frame(read_csv("latent_factors_M3_wide.csv"))
factor_mofa <- factor_mofa %>% select(sample,Factor7)
factor_mofa$sample <- dataset$SampleID
factor_mofa <- dataset_met %>% left_join(factor_mofa,by=c("SampleID"="sample"))
factor_metab_corr <- partial_spearman_corr(x_df = factor_mofa,y_df = dataset_met,covar_indices = 2:7,x_indices = 88,y_indices = 8:87,n_cores = 8,thresh = 0.2,include_within = FALSE)
factor_mgs_corr <- partial_spearman_corr(x_df = factor_mofa,y_df = combined$MGS,covar_indices = 2:7,x_indices = 88,y_indices = 2:210,n_cores = 8,thresh = 0.2,include_within = FALSE)
