library(readr)
library(tidyverse)
library(VIM)
library(caret)
library(qvalue)
dataset <- as.data.frame(read_csv("combined_df_M3.csv",show_col_types = FALSE))
metagenome <- as.data.frame(read_csv("metagenome combined.csv"))
metabolites <- as.data.frame(read_csv("combined_df_M1.csv"))
metabolites <- metabolites %>% select(1:7,14:176)
bugs <- as.data.frame(read_csv("Cluster 12-13- Metabolite Associated.csv"))
MGS_interest <- metagenome %>% select(2,bugs$name)
dataset_MGS <- dataset %>% select(1,2:8) %>% left_join(MGS_interest,by=c("SampleID"="MGSampleID"))
run_lm_parallel <- function(df, x_indices, y_indices, covar_indices, seed = 42, n_cores = 8) {
  library(parallel)
  library(dplyr)
  library(qvalue)  # Load qvalue package
  
  set.seed(seed)
  
  combos <- expand.grid(x = x_indices, y = y_indices)
  
  results <- mclapply(1:nrow(combos), function(i) {
    x_idx <- combos$x[i]
    y_idx <- combos$y[i]
    
    x_name <- colnames(df)[x_idx]
    y_name <- colnames(df)[y_idx]
    
    model_data <- df[, c(x_idx, y_idx, covar_indices)]
    colnames(model_data)[1:2] <- c("x", "y")
    model_data <- na.omit(model_data)
    
    if (nrow(model_data) < length(covar_indices) + 3) return(NULL)
    
    formula <- as.formula(paste("y ~ x +", paste(colnames(df)[covar_indices], collapse = " + ")))
    
    fit <- tryCatch({
      lm(formula, data = model_data)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(fit)) {
      summary_fit <- summary(fit)
      p_val <- coef(summary_fit)["x", "Pr(>|t|)"]
      estimate <- coef(summary_fit)["x", "Estimate"]
      
      return(data.frame(
        x_var = x_name,
        y_var = y_name,
        estimate = estimate,
        p_value = p_val
      ))
    } else {
      return(NULL)
    }
  }, mc.cores = n_cores)
  
  results_df <- bind_rows(results)
  
  if (nrow(results_df) > 1) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = "fdr")
    
    # Compute q-values
    qobj <- qvalue(p = results_df$p_value)
    results_df$q_value <- qobj$qvalues
  } else {
    results_df$p_adj <- results_df$p_value
    results_df$q_value <- results_df$p_value
  }
  
  return(results_df)
}
run_partial_spearman_parallel <- function(df, x_indices, y_indices, covar_indices, seed = 42, n_cores = 8) {
  library(parallel)
  library(dplyr)
  library(ppcor)
  library(qvalue)
  
  set.seed(seed)
  
  combos <- expand.grid(x = x_indices, y = y_indices)
  
  results <- mclapply(1:nrow(combos), function(i) {
    x_idx <- combos$x[i]
    y_idx <- combos$y[i]
    
    x_name <- colnames(df)[x_idx]
    y_name <- colnames(df)[y_idx]
    
    data_subset <- df[, c(x_idx, y_idx, covar_indices)]
    colnames(data_subset)[1:2] <- c("x", "y")
    data_subset <- na.omit(data_subset)
    
    if (nrow(data_subset) < length(covar_indices) + 3) return(NULL)
    
    test_result <- tryCatch({
      pcor.test(data_subset$x, data_subset$y, data_subset[, -c(1, 2)], method = "spearman")
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(test_result)) {
      return(data.frame(
        x_var = x_name,
        y_var = y_name,
        estimate = test_result$estimate,
        p_value = test_result$p.value
      ))
    } else {
      return(NULL)
    }
  }, mc.cores = n_cores)
  
  results_df <- bind_rows(results)
  
  if (nrow(results_df) > 1) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = "fdr")
    qobj <- qvalue(p = results_df$p_value)
    results_df$q_value <- qobj$qvalues
  } else {
    results_df$p_adj <- results_df$p_value
    results_df$q_value <- results_df$p_value
  }
  
  return(results_df)
}
lm_metformin_MGS <- run_lm_parallel(dataset_MGS,x_indices = 8,y_indices = 9:114,covar_indices = 2:7)
corr_metformin_MGS <- run_partial_spearman_parallel(dataset_MGS,x_indices = 8,y_indices = 9:114,covar_indices = 2:7)
write.csv(lm_metformin_MGS,"lm(n) met MGS.csv")
write.csv(corr_metformin_MGS,"corr met MGS.csv")
corr_metformin_MGS_sig <- corr_metformin_MGS %>% filter(p_adj < 0.05)
dataset_MGS_met <- metagenome %>% select(2, corr_metformin_MGS_sig$y_var) %>% rename(SampleID = MGSampleID)
metagenome_other <- as.data.frame(read_csv("metagenome without cags.csv"))
run_lm_by_datatype <- function(x_df, y_df, covar_df,
                               kegg_ko_indices,
                               kegg_module_indices,
                               kegg_pathway_indices,
                               gmm_indices,
                               x_combinations = NULL) {
  library(dplyr)
  library(parallel)
  library(purrr)
  
  # Dynamically detect first column names
  id_x <- colnames(x_df)[1]
  id_y <- colnames(y_df)[1]
  id_covar <- colnames(covar_df)[1]
  
  # Rename to a common join key
  x_df <- x_df %>% rename(JoinID = !!id_x)
  y_df <- y_df %>% rename(JoinID = !!id_y)
  covar_df <- covar_df %>% rename(JoinID = !!id_covar)
  
  # Merge all data
  merged_df <- x_df %>%
    inner_join(y_df, by = "JoinID") %>%
    inner_join(covar_df, by = "JoinID")
  
  x_vars <- setdiff(colnames(x_df), "JoinID")
  covar_vars <- setdiff(colnames(covar_df), "JoinID")
  
  data_types <- list(
    KEGG_KO = kegg_ko_indices,
    KEGG_Module = kegg_module_indices,
    KEGG_Pathway = kegg_pathway_indices,
    GMM = gmm_indices
  )
  
  convert_x_combinations <- function(combos, x_vars) {
    lapply(combos, function(group) {
      if (is.numeric(group)) {
        return(x_vars[group])
      } else if (all(group %in% x_vars)) {
        return(group)
      } else {
        stop("Invalid x_combinations: must be indices or valid x_var names.")
      }
    })
  }
  
  run_models_for_data_type <- function(data_type, y_indices) {
    y_vars <- colnames(y_df)[y_indices]
    results <- list()
    
    for (x_var in x_vars) {
      for (y_var in y_vars) {
        formula <- as.formula(
          paste(y_var, "~", x_var, "+", paste(covar_vars, collapse = " + "))
        )
        
        fit <- tryCatch({
          lm(formula, data = merged_df)
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          summary_fit <- summary(fit)
          if (x_var %in% rownames(coef(summary_fit))) {
            p_val <- coef(summary_fit)[x_var, "Pr(>|t|)"]
            estimate <- coef(summary_fit)[x_var, "Estimate"]
            
            results[[length(results) + 1]] <- data.frame(
              x_var = x_var,
              y_var = y_var,
              data_type = data_type,
              estimate = estimate,
              p_value = p_val
            )
          }
        }
      }
    }
    
    if (length(results) > 0) {
      df <- dplyr::bind_rows(results)
      df$p_adj <- p.adjust(df$p_value, method = "fdr")
      return(df)
    } else {
      return(NULL)
    }
  }
  
  os_type <- .Platform$OS.type
  num_cores <- detectCores() - 1
  
  if (os_type == "windows") {
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist = c("x_vars", "covar_vars", "merged_df", "y_df", "run_models_for_data_type"), envir = environment())
    all_results <- parLapply(cl, names(data_types), function(dt) {
      run_models_for_data_type(dt, data_types[[dt]])
    })
    stopCluster(cl)
  } else {
    all_results <- mclapply(names(data_types), function(dt) {
      run_models_for_data_type(dt, data_types[[dt]])
    }, mc.cores = num_cores)
  }
  
  names(all_results) <- names(data_types)
  all_results <- all_results[!sapply(all_results, is.null)]
  
  combined_df <- dplyr::bind_rows(all_results)
  significant_df <- combined_df %>% filter(p_adj < 0.05)
  insignificant_df <- combined_df %>% filter(p_adj >= 0.05)
  
  summary_df <- combined_df %>%
    group_by(x_var, data_type) %>%
    summarise(
      total = n(),
      significant = sum(p_adj < 0.05),
      percent_significant = round(100 * significant / total, 1),
      .groups = "drop"
    )
  
  compute_custom_combinations_union <- function(sig_df, x_combinations, x_vars) {
    x_combinations <- convert_x_combinations(x_combinations, x_vars)
    data_types <- unique(sig_df$data_type)
    y_counts <- sig_df %>%
      group_by(data_type) %>%
      summarise(total_y = n_distinct(y_var), .groups = "drop")
    
    map_dfr(x_combinations, function(x_group) {
      map_dfr(data_types, function(dt) {
        y_union <- sig_df %>%
          filter(x_var %in% x_group, data_type == dt) %>%
          pull(y_var) %>%
          unique()
        
        total_y <- y_counts %>%
          filter(data_type == dt) %>%
          pull(total_y)
        
        tibble(
          x_vars = paste(x_group, collapse = " + "),
          data_type = dt,
          significant_union = length(y_union),
          total_y = total_y,
          proportion_significant = round(length(y_union) / total_y, 3)
        )
      })
    })
  }
  
  custom_union_df <- NULL
  if (!is.null(x_combinations)) {
    custom_union_df <- compute_custom_combinations_union(significant_df, x_combinations, x_vars)
  }
  
  return(list(
    summary = summary_df,
    significant_results = significant_df,
    insignificant_results = insignificant_df,
    custom_union = custom_union_df
  ))
}
covar_df <- dataset %>% select(1:7)
associated_features <- as.data.frame(read_csv("combined_df_M1.csv"))
associated_features_colnames <- colnames(associated_features)[12:848]
associated_features_colnames <- intersect(colnames(metagenome_other),associated_features_colnames)
metagenome_other <- metagenome_other %>% select(2,all_of(associated_features_colnames))
metagenome_other <- metagenome_other %>% rename('SampleID'='MGSampleID')
associated_microbes_kegg_datasets <- run_lm_by_datatype(kegg_ko_indices = 2:440,kegg_pathway_indices = 441:461,kegg_module_indices = 462:514,gmm_indices = 515:516,x_df = dataset_MGS_met,covar_df = covar_df,y_df = metagenome_other,x_combinations =list(c(2,3,10),c(13,14)))
annotation <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.variables.v6.r.csv"))
annotation <- annotation %>% select(VariableID,DisplayName,Notes)
extract_phylum_and_microbe <- function(df) {
  df$Phylum <- sapply(df$Notes, function(tax_string) {
    parts <- unlist(strsplit(tax_string, ";"))
    if (length(parts) >= 2) {
      return(parts[2])
    } else {
      return(NA)
    }
  })
  
  df$Microbe <- sapply(df$Notes, function(tax_string) {
    parts <- unlist(strsplit(tax_string, ";"))
    if (length(parts) >= 7) {
      return(parts[7])  # Assuming the 7th level is the species or strain
    } else if (length(parts) >= 6) {
      return(parts[6])  # Fallback to genus if species is missing
    } else {
      return(NA)
    }
  })
  
  return(df)
}
# Apply the function
annotation <- extract_phylum_and_microbe(annotation)
corr_metformin_MGS_sig <- corr_metformin_MGS_sig %>%
  left_join(annotation,by=c("y_var"="VariableID"))
corr_metformin_MGS_sig_out <- corr_metformin_MGS_sig %>% select(estimate,p_value,p_adj,DisplayName,Phylum,Microbe.x)
corr_metformin_MGS <- left_join(corr_metformin_MGS,annotation,by=c("y_var"="VariableID"))
corr_metformin_MGS_out <- corr_metformin_MGS %>% select(estimate,p_value,p_adj,DisplayName.x,Phylum.x,Microbe)
write_csv(corr_metformin_MGS_sig_out,"corr_annotated_sig.csv")
write_csv(corr_metformin_MGS_out,"corr_annotated.csv")
