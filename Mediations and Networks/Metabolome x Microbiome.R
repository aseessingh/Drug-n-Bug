library(readr)
library(dplyr)
metabolite_df <- as.data.frame(read_csv("metabolites.csv"))
metabolite_df <- metabolite_df %>% select(-c(1,3:8))
microbe_df <-  as.data.frame(read_csv("MGS.csv"))
microbe_df <- microbe_df %>% select(-1)
covar_df <- as.data.frame(read_csv("covar df.csv"))
covar_df <- covar_df %>% select(-1)
corr_df <- as.data.frame(read_csv("corr met MGS.csv",show_col_types = FALSE))
bugs_of_interest <- corr_df %>% filter(p_adj<0.05) %>% select(y_var)
run_simple_mediation_pipeline <- function(
    metabolite_df, microbe_df, covar_df, bugs_of_interest_df, graphml_path,
    n_boot = 1000, n_cores = 8, output_prefix = "simple_output"
) {
  library(glmnet)
  library(igraph)
  library(tidyverse)
  library(parallel)
  
  # Rename first column to 'sample_id'
  metabolite_df <- metabolite_df %>% rename(sample_id = 1)
  microbe_df <- microbe_df %>% rename(sample_id = 1)
  covar_df <- covar_df %>% rename(sample_id = 1)
  
  # Merge all data frames
  data <- reduce(list(metabolite_df, microbe_df, covar_df), full_join, by = "sample_id")
  rownames(data) <- data$sample_id
  data <- data %>% select(-sample_id)
  
  metabolites <- metabolite_df %>% select(-sample_id)
  microbes <- microbe_df %>% select(-sample_id)
  covars <- covar_df %>% select(-sample_id)
  covar_names <- colnames(covars)
  bugs_of_interest <- bugs_of_interest_df[[1]]
  
  # Standardize microbes and covariates
  microbes <- as.data.frame(scale(microbes))
  covars <- as.data.frame(scale(covars))
  
  # Linear regression
  linear_results <- list()
  for (metab in colnames(metabolites)) {
    y <- metabolites[[metab]]
    for (microbe in colnames(microbes)) {
      x <- microbes[[microbe]]
      fit <- lm(y ~ x + ., data = covars)
      summary_fit <- summary(fit)
      beta <- summary_fit$coefficients[2, 1]
      p <- summary_fit$coefficients[2, 4]
      linear_results[[length(linear_results) + 1]] <- data.frame(
        metabolite = metab,
        microbe = microbe,
        beta = beta,
        pvalue = p
      )
    }
  }
  linear_df <- bind_rows(linear_results)
  linear_df$p_adj <- p.adjust(linear_df$pvalue, method = "fdr")
  significant_linear <- linear_df %>% filter(p_adj < 0.05)
  write_csv(significant_linear, paste0(output_prefix, "_significant_linear.csv"))
  
  # LASSO regression with bootstrapping
  bootstrap_lasso <- function(metab_name) {
    y <- metabolites[[metab_name]]
    
    run_bootstrap <- function(i) {
      idx <- sample(seq_along(y), replace = TRUE)
      x_boot <- as.matrix(cbind(microbes[idx, ], covars[idx, ]))
      y_boot <- y[idx]
      cvfit <- cv.glmnet(x_boot, y_boot, alpha = 1)
      model <- glmnet(x_boot, y_boot, alpha = 1, lambda = cvfit$lambda.min)
      coefs <- coef(model)
      nonzero <- which(coefs != 0)[-1]  # exclude intercept
      if (length(nonzero) == 0) return(NULL)
      data.frame(
        microbe = rownames(coefs)[nonzero],
        beta = as.numeric(coefs[nonzero])
      )
    }
    
    boot_results <- mclapply(
      1:n_boot, run_bootstrap,
      mc.cores = min(8, detectCores()),
      mc.preschedule = FALSE
    )
    boot_results <- boot_results[!sapply(boot_results, is.null)]
    all_coefs <- bind_rows(boot_results)
    
    if (nrow(all_coefs) == 0) return(NULL)
    
    summary_df <- all_coefs %>%
      group_by(microbe) %>%
      summarise(
        freq = n() / n_boot,
        avg_beta = mean(beta),
        .groups = "drop"
      ) %>%
      mutate(metabolite = metab_name)
    
    summary_df
  }
  
  selected_microbes <- mclapply(
    colnames(metabolites),
    bootstrap_lasso,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  selected_microbes <- selected_microbes[!sapply(selected_microbes, is.null)]
  lasso_df <- bind_rows(selected_microbes)
  
  # Save all selection frequencies
  write_csv(lasso_df, paste0(output_prefix, "_lasso_frequencies.csv"))
  
  # Filter significant ones
  lasso_df$p_adj <- p.adjust(1 - lasso_df$freq, method = "fdr")
  significant_lasso <- lasso_df %>% filter(p_adj < 0.05)
  write_csv(significant_lasso, paste0(output_prefix, "_significant_lasso.csv"))
  
  # Network filtering
  g <- read_graph(graphml_path, format = "graphml")
  all_significant_microbes <- unique(c(significant_linear$microbe, significant_lasso$microbe))
  
  linked_microbes <- all_significant_microbes[
    sapply(all_significant_microbes, function(microbe) {
      any(sapply(bugs_of_interest, function(bug) {
        bug %in% V(g)$name && microbe %in% V(g)$name &&
          length(all_simple_paths(g, from = bug, to = microbe)) > 0
      }))
    })
  ]
  
  write_csv(data.frame(microbe = linked_microbes), paste0(output_prefix, "_linked_microbes.csv"))
  
  return(list(
    significant_linear = significant_linear,
    significant_lasso = significant_lasso,
    lasso_frequencies = lasso_df,
    linked_microbes = linked_microbes
  ))
}
metabxmicrobe <- run_simple_mediation_pipeline(metabolite_df,microbe_df,covar_df,graphml_path = "correlation_network.graphml",bugs_of_interest,n_cores = 128)
combine_linear_lasso_with_links <- function(
    lasso_csv_path,
    linear_csv_path,
    graphml_path,
    bugs_of_interest_df,
    covar_df,
    output_csv_path = "linked_microbes_with_sources.csv",
    freq_cutoff = 0.75
) {
  library(dplyr)
  library(readr)
  library(igraph)
  library(purrr)
  
  # Extract covariate names
  covariate_names <- colnames(covar_df)[2:7]
  
  # Load linear results and extract microbes
  linear_df <- read_csv(linear_csv_path, show_col_types = FALSE)
  linear_microbes <- linear_df$microbe
  
  # Load LASSO results and filter
  lasso_df <- read_csv(lasso_csv_path, show_col_types = FALSE)
  filtered_lasso <- lasso_df %>%
    filter(freq >= freq_cutoff & !(microbe %in% covariate_names))
  
  # Combine unique microbes
  combined_microbes <- unique(c(filtered_lasso$microbe, linear_microbes))
  combined_df <- filtered_lasso %>% filter(microbe %in% combined_microbes)
  
  # Load graph and clean inputs
  g <- read_graph(graphml_path, format = "graphml")
  graph_nodes <- V(g)$name
  bugs <- intersect(bugs_of_interest_df[[1]], graph_nodes)
  candidates <- intersect(combined_df$microbe, graph_nodes)
  
  # Compute connectivity: for each microbe, list bugs it's linked to
  linked_info <- map_df(candidates, function(microbe) {
    connected_bugs <- bugs[sapply(bugs, function(bug) {
      distances(g, v = bug, to = microbe)[1, 1] < Inf
    })]
    if (length(connected_bugs) == 0) return(NULL)
    tibble(
      microbe = microbe,
      connected_to = paste(connected_bugs, collapse = "; ")
    )
  })
  
  # Final output
  final_df <- combined_df %>%
    inner_join(linked_info, by = "microbe")
  
  write_csv(final_df, output_csv_path)
  return(final_df)
}
result <- combine_linear_lasso_with_links(
  lasso_csv_path = "simple_output_lasso_frequencies.csv",
  linear_csv_path = "simple_output_significant_linear.csv",
  graphml_path = "correlation_network.graphml",
  bugs_of_interest_df = bugs_of_interest,
  covar_df = covar_df,
  output_csv_path = "linked_microbes_combined.csv"
)


