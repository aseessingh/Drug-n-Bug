run_glmnet_model_pipeline <- function(
    data_file,
    mediator_file,
    hba1c_index,
    ag_index,
    dose_index,
    serum_index,
    covariate_indices,
    output_csv = "glmnet_model_summary.csv",
    model_rds = "best_glmnet_model.rds",
    plot_file = "model_performance_plot.png",
    n_boot = 1000,
    n_cores = parallel::detectCores() - 1
) {
  library(glmnet)
  library(caret)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(readr)
  library(ggplot2)
  
  set.seed(123)
  
  # Read input data
  df <- as.data.frame(read_csv(data_file))
  mediator_df <- read_csv(mediator_file)
  
  # Extract mediator variable names
  mediator_names <- mediator_df$mediator
  mediator_indices <- which(names(df) %in% mediator_names)
  
  # Extract variable names
  target_var <- names(df)[dose_index]
  feature_indices <- unique(c(hba1c_index, ag_index, serum_index, covariate_indices, mediator_indices))
  feature_names <- names(df)[feature_indices]
  
  # Subset and prepare data
  df_model <- df[, c(target_var, feature_names)]
  df_model <- na.omit(df_model)
  
  # Split data
  trainIndex <- createDataPartition(df_model[[target_var]], p = 0.7, list = FALSE)
  train <- df_model[trainIndex, ]
  temp <- df_model[-trainIndex, ]
  valIndex <- createDataPartition(temp[[target_var]], p = 0.5, list = FALSE)
  validation <- temp[valIndex, ]
  test <- temp[-valIndex, ]
  
  # Model matrices
  formula <- as.formula(paste(target_var, "~ ."))
  x_train <- model.matrix(formula, train)[, -1]
  y_train <- train[[target_var]]
  x_val <- model.matrix(formula, validation)[, -1]
  y_val <- validation[[target_var]]
  x_test <- model.matrix(formula, test)[, -1]
  y_test <- test[[target_var]]
  
  # Define finer alpha and lambda grids
  alphas <- seq(0, 1, by = 0.005)
  lambdas <- 10^seq(3, -4, length = 200)
  
  # Bootstrap function
  bootstrap_model <- function(i) {
    a <- alphas[(i %% length(alphas)) + 1]
    boot_idx <- sample(1:nrow(x_train), replace = TRUE)
    x_boot <- x_train[boot_idx, ]
    y_boot <- y_train[boot_idx]
    cv_model <- cv.glmnet(x_boot, y_boot, alpha = a, lambda = lambdas, nfolds = 10)
    
    pred_val <- predict(cv_model, newx = x_val, s = "lambda.min")
    rmsee <- sqrt(mean((pred_val - y_val)^2))
    mae <- mean(abs(pred_val - y_val))
    r2 <- 1 - sum((pred_val - y_val)^2) / sum((y_val - mean(y_val))^2)
    
    dropped <- colnames(x_train)[which(coef(cv_model, s = "lambda.min") == 0)]
    list(model = cv_model, alpha = a, lambda = cv_model$lambda.min,
         rmsee = rmsee, mae = mae, r2 = r2,
         dropped = dropped[!dropped %in% target_var])
  }
  
  # Run in parallel
  results <- mclapply(1:n_boot, bootstrap_model, mc.cores = n_cores)
  
  # Compile results
  df_results <- bind_rows(lapply(seq_along(results), function(i) {
    res <- results[[i]]
    data.frame(model_id = i, alpha = res$alpha, lambda = res$lambda,
               rmsee = res$rmsee, mae = res$mae, r2 = res$r2,
               dropped = paste(res$dropped, collapse = ";"), stringsAsFactors = FALSE)
  }))
  
  # Normalize metrics and compute composite score
  df_results <- df_results %>%
    mutate(
      norm_rmsee = (rmsee - min(rmsee)) / (max(rmsee) - min(rmsee)),
      norm_mae = (mae - min(mae)) / (max(mae) - min(mae)),
      norm_r2 = (r2 - min(r2)) / (max(r2) - min(r2)),
      composite_score = norm_r2 - norm_rmsee - norm_mae
    ) %>%
    arrange(desc(composite_score)) %>%
    slice(1:50)
  
  # Evaluate top 50 on test set
  final_results <- df_results %>%
    rowwise() %>%
    mutate(test_pred = list(predict(results[[model_id]]$model, newx = x_test, s = results[[model_id]]$lambda)),
           test_rmsee = sqrt(mean((test_pred[[1]] - y_test)^2)),
           test_mae = mean(abs(test_pred[[1]] - y_test)),
           test_r2 = 1 - sum((test_pred[[1]] - y_test)^2) / sum((y_test - mean(y_test))^2)) %>%
    select(model_id, alpha, lambda, dropped, test_rmsee, test_mae, test_r2)
  
  # Save output summary
  write.csv(final_results, output_csv, row.names = FALSE)
  
  # Save best model
  best_model_id <- df_results$model_id[1]
  best_model <- results[[best_model_id]]$model
  saveRDS(best_model, file = model_rds)
  
  # Plotting
  best_model_row <- final_results[1, ]
  p <- ggplot(final_results, aes(x = test_rmsee, y = test_mae, color = test_r2)) +
    geom_point(size = 3) +
    scale_color_viridis_c(name = "Test R²") +
    geom_point(data = best_model_row, aes(x = test_rmsee, y = test_mae), color = "red", size = 4, shape = 21, fill = "white") +
    geom_text(data = best_model_row, aes(label = "Best Model"), vjust = -1, color = "red") +
    labs(title = "Model Performance: RMSEE vs MAE (colored by R²)",
         x = "Test RMSEE", y = "Test MAE") +
    theme_minimal()
  
  ggsave(plot_file, plot = p, width = 10, height = 6)
}
run_glmnet_model_pipeline(
  data_file = "combined_df_M1.csv",
  mediator_file = "sig_M1.csv",
  hba1c_index = 10,
  ag_index = 11,
  dose_index = 8,
  serum_index = 9,
  covariate_indices = 2:7,
  output_csv = "glmnet_M1.csv",
  model_rds = "glmnet_M1_model.rds",
  plot_file = "m1_perf.png"
)
run_glmnet_model_pipeline(
  data_file = "combined_df_M3.csv",
  mediator_file = "sig_M3.csv",
  hba1c_index = 10,
  ag_index = 11,
  dose_index = 8,
  serum_index = 9,
  covariate_indices = 2:7,
  output_csv = "glmnet_M3.csv",
  model_rds = "glmnet_M3_model.rds",
  plot_file = "m3_perf.png"
)

