library(readr)
args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
total_chunks <- as.integer(args[2])

M1 <- as.data.frame(read_csv("combined_df_M1.csv"))
M3 <- as.data.frame(read_csv("combined_df_M3.csv"))

#Mediation
run_extended_mediation <- function(data, input_cols, output_cols, mediator_cols, covar_cols = NULL,
                                   n_boot = 1000, seed = 123, num_cores = 8, x = "results",
                                   output_dir = "Chunked Results") {
  library(mediation)
  library(parallel)
  library(dplyr)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  RNGkind("L'Ecuyer-CMRG")  # For reproducible parallel RNG
  set.seed(seed)
  
  input_cols <- colnames(data)[input_cols]
  output_cols <- colnames(data)[output_cols]
  mediator_cols <- colnames(data)[mediator_cols]
  covar_cols <- if (!is.null(covar_cols)) colnames(data)[covar_cols] else character(0)
  
  combos <- expand.grid(input_cols, output_cols, stringsAsFactors = FALSE)
  
  model_summaries <- list()  # To store all model summaries and diagnostics
  
  results <- mclapply(mediator_cols, function(m) {
    causal_results <- list()
    
    for (i in 1:nrow(combos)) {
      x_var <- combos[i, 1]
      y_var <- combos[i, 2]
      
      required_cols <- c(x_var, y_var, m, covar_cols)
      if (!all(required_cols %in% colnames(data))) {
        message(paste("Skipping due to missing columns:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
        next
      }
      
      df <- data[, required_cols]
      colnames(df)[1:3] <- c("x", "y", "m")
      
      try({
        message(paste("Processing:", m, "->", x_var, "→", y_var))
        
        # Causal Mediation
        med_model <- lm(m ~ x + ., data = df)
        out_model <- lm(y ~ m + x + ., data = df)
        
        # Extract and store model summaries and diagnostics
        extract_model_info <- function(model, model_type) {
          s <- summary(model)
          df <- as.data.frame(coef(s))
          df$model <- model_type
          df$input <- x_var
          df$output <- y_var
          df$mediator <- m
          df$term <- rownames(df)
          df$r.squared <- s$r.squared
          df$adj.r.squared <- s$adj.r.squared
          df$residual.se <- s$sigma
          df$f.statistic <- s$fstatistic[1]
          df$f.df1 <- s$fstatistic[2]
          df$f.df2 <- s$fstatistic[3]
          df
        }
        
        med_df <- extract_model_info(med_model, "mediator")
        out_df <- extract_model_info(out_model, "outcome")
        
        model_summaries[[length(model_summaries) + 1]] <<- bind_rows(med_df, out_df)
        
        # Mediation analysis
        med_out <- tryCatch({
          mediate(med_model, out_model, treat = "x", mediator = "m", boot = TRUE, sims = n_boot)
        }, error = function(e) {
          message(paste("Causal mediation failed for", m, "->", x_var, "→", y_var, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(med_out)) {
          causal_results[[length(causal_results) + 1]] <- data.frame(
            input = x_var, output = y_var, mediator = m,
            acme = med_out$d0, acme_ci_lower = med_out$d0.ci[1], acme_ci_upper = med_out$d0.ci[2],
            ade = med_out$z0, ade_ci_lower = med_out$z0.ci[1], ade_ci_upper = med_out$z0.ci[2],
            total_effect = med_out$tau.coef, total_ci_lower = med_out$tau.ci[1], total_ci_upper = med_out$tau.ci[2],
            prop_mediated = med_out$n0,
            acme_p = med_out$d0.p, ade_p = med_out$z0.p, total_p = med_out$tau.p
          )
        }
      }, silent = TRUE)
    }
    
    list(causal = causal_results)
  }, mc.cores = num_cores)
  
  # Combine results
  causal_df <- bind_rows(lapply(results, function(r) bind_rows(r$causal)))
  
  # Adjust p-values
  if ("acme_p" %in% colnames(causal_df)) {
    causal_df <- causal_df %>%
      group_by(input, output) %>%
      mutate(acme_p_adj = p.adjust(acme_p, method = "BH"),
             ade_p_adj = p.adjust(ade_p, method = "BH"),
             total_p_adj = p.adjust(total_p, method = "BH")) %>%
      ungroup()
  }
  
  # Write causal results
  write.csv(causal_df, file.path(output_dir, paste0(x, "_causal.csv")), row.names = FALSE)
  
  # Write model summaries and diagnostics
  if (length(model_summaries) > 0) {
    summary_df <- bind_rows(model_summaries)
    write.csv(summary_df, file.path(output_dir, paste0(x, "_model_summaries.csv")), row.names = FALSE)
  }
}

get_chunk_indices <- function(all_mediators, chunk_id, total_chunks) {
  chunk_size <- ceiling(length(all_mediators) / total_chunks)
  start_idx <- (chunk_id - 1) * chunk_size + 1
  end_idx <- min(chunk_id * chunk_size, length(all_mediators))
  return(all_mediators[start_idx:end_idx])
}
mediators_M1 <- 12:848
mediators_M3 <- 12:221
mediator_chunk_M1 <- get_chunk_indices(mediators_M1, chunk_id, total_chunks)
mediator_chunk_M3 <- get_chunk_indices(mediators_M3, chunk_id, total_chunks)
run_extended_mediation(
  data = M1,
  input_cols = 8:9,
  output_cols = 10:11,
  mediator_cols = mediator_chunk_M1,
  num_cores = 8,
  covar_cols = 2:7,
  x = paste0("mediation_M1_chunk_", chunk_id)
)
run_extended_mediation(
  data = M3,
  input_cols = 8:9,
  output_cols = 10:11,
  mediator_cols = mediator_chunk_M3,
  num_cores = 8,
  covar_cols = 2:7,
  x = paste0("mediation_M3_chunk_", chunk_id)
)


