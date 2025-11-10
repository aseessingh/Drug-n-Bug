library(stats)
library(dplyr)
library(readr)
library(parallel)

# Input
metabolite <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv"))
metadata <- read.csv("metadata3.csv")
results_metab <- read.csv("drug_non_drug_base_metab_stats.csv")
results_mgs <- read.csv("drug_non_drug_base_mgs_stats.csv")
results_kegg1 <- read.csv("drug_non_drug_base_KEGG_KO_stats.csv")
results_kegg2 <- read.csv("drug_non_drug_base_KEGG_module_stats.csv")
results_gmm <- read.csv("drug_non_drug_base_GMM_stats.csv")
results_kegg3 <- read.csv("drug_non_drug_base_KEGG_pathway_stats.csv")
mgs <- read.csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv")
gmm <- read.csv("hub.function_adjusted.GMM.down.10000000.fpkm.mgsamples.v5.r.csv")
kegg_1 <- read.csv("hub.function_adjusted.KEGG_ko.down.10000000.fpkm.mgsamples.v3.r.csv")
kegg_2 <- read.csv("hub.function_adjusted.KEGG_module.down.10000000.fpkm.mgsamples.v3.r.csv")
kegg_3 <- read.csv("hub.function_adjusted.KEGG_pathway.down.10000000.fpkm.mgsamples.v3.r.csv")
group3 <- read.csv("group 3.csv")

# Preprocessing
group3 <- group3 %>% filter(SampleID %in% metadata$SampleID) %>% dplyr::select(SampleID, GLYCATHB)
metadata <- left_join(metadata, group3, by = "SampleID")
metabolite_AG <- metabolite %>% dplyr::select(MGS_PID, `1,5-anhydroglucitol (1,5-AG)`) %>% rename(AG = `1,5-anhydroglucitol (1,5-AG)`)
metadata <- left_join(metadata, metabolite_AG, by = c("SampleID" = "MGS_PID"))
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
metabolite_results <- metabolite %>% dplyr::select(MGS_PID, results_metab$Feature)
gmm_results <- gmm %>% dplyr::select(MGSampleID,results_gmm$Feature)
kegg1_results <- kegg_1 %>% dplyr::select(MGSampleID,results_kegg1$Feature)
kegg2_results <- kegg_2 %>% dplyr::select(MGSampleID,results_kegg2$Feature)
kegg3_results <- kegg_3 %>% dplyr::select(MGSampleID,results_kegg3$Feature)
mgs_results <- mgs %>% dplyr::select(MGSampleID,results_mgs$Feature)
combined_df_metab <- left_join(metadata, metabolite_results, by = c("SampleID" = "MGS_PID"))
combined_df_gmm <- left_join(metadata,gmm_results,by=c("SampleID"="MGSampleID"))
combined_df_kegg <- left_join(metadata,kegg1_results,by=c("SampleID"="MGSampleID"))
combined_df_kegg <- left_join(combined_df_kegg,kegg2_results,by=c("SampleID"="MGSampleID"))
combined_df_kegg <- left_join(combined_df_kegg,kegg3_results,by=c("SampleID"="MGSampleID"))
combined_df_mgs <- left_join(metadata,mgs_results,by=c("SampleID"="MGSampleID"))
# Mediation
run_mediation_analysis <- function(data, input_cols, output_cols, mediator_cols, covar_cols, n_boot = 1000, seed = 123, numOfCores = 1) {
  set.seed(seed)
  
  sanitize_name <- function(name) {
    gsub("[^[:alnum:]_]", "_", name)
  }
  
  results <- list()
  covar_names <- colnames(data)[covar_cols]
  
  for (x_col in input_cols) {
    for (y_col in output_cols) {
      for (m_col in mediator_cols) {
        x_name <- colnames(data)[x_col]
        y_name <- colnames(data)[y_col]
        m_name <- colnames(data)[m_col]
        mediator_label <- sanitize_name(m_name)
        
        required_cols <- c(x_name, y_name, m_name, covar_names)
        
        if (!all(required_cols %in% colnames(data))) {
          warning("Skipping due to missing columns: ", paste(setdiff(required_cols, colnames(data)), collapse = ", "))
          next
        }
        
        df <- data[, required_cols]
        colnames(df)[1:3] <- c("x", "y", "m")
        
        model_m <- tryCatch(lm(m ~ x + ., data = df), error = function(e) NULL)
        model_y <- tryCatch(lm(y ~ x + m + ., data = df), error = function(e) NULL)
        model_total <- tryCatch(lm(y ~ x + ., data = df), error = function(e) NULL)
        
        if (is.null(model_m) || is.null(model_y) || is.null(model_total)) next
        if (!"m" %in% rownames(coef(summary(model_y)))) {
          warning(paste("Mediator 'm' not found in model_y for", m_name, "- skipping."))
          next
        }
        
        a <- coef(summary(model_m))["x", "Estimate"]
        b <- coef(summary(model_y))["m", "Estimate"]
        ab <- a * b
        ade <- coef(summary(model_y))["x", "Estimate"]
        c_total <- coef(summary(model_total))["x", "Estimate"]
        prop_med <- ifelse(c_total != 0, ab / c_total, NA)
        
        se_a <- coef(summary(model_m))["x", "Std. Error"]
        se_b <- coef(summary(model_y))["m", "Std. Error"]
        se_ab <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)
        z <- ab / se_ab
        p_sobel <- 2 * (1 - pnorm(abs(z)))
        
        ab_boot <- parallel::mclapply(1:n_boot, function(i) {
          idx <- sample(seq_len(nrow(df)), replace = TRUE)
          boot_df <- df[idx, ]
          boot_model_m <- tryCatch(lm(m ~ x + ., data = boot_df), error = function(e) NULL)
          boot_model_y <- tryCatch(lm(y ~ x + m + ., data = boot_df), error = function(e) NULL)
          if (is.null(boot_model_m) || is.null(boot_model_y)) return(NA)
          if (!all(c("x", "m") %in% names(coef(boot_model_m))) || !"m" %in% names(coef(boot_model_y))) return(NA)
          a_b <- coef(boot_model_m)["x"]
          b_b <- coef(boot_model_y)["m"]
          a_b * b_b
        }, mc.cores = numOfCores)
        
        ab_boot <- unlist(ab_boot)
        ab_boot <- ab_boot[!is.na(ab_boot)]
        
        ci_lower <- quantile(ab_boot, 0.025, na.rm = TRUE)
        ci_upper <- quantile(ab_boot, 0.975, na.rm = TRUE)
        p_boot <- mean(ab_boot < 0 | ab_boot > 0) * 2
        
        results[[length(results) + 1]] <- data.frame(
          input = x_name,
          output = y_name,
          mediator = mediator_label,
          a = a,
          b = b,
          acme = ab,
          ade = ade,
          total_effect = c_total,
          prop_mediated = prop_med,
          p_sobel = p_sobel,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          p_boot = p_boot
        )
      }
    }
  }
  
  result_df <- dplyr::bind_rows(results)
  result_df$p_adj_boot <- p.adjust(result_df$p_boot, method = "BH")
  return(result_df)
}
# Process all mediators
mediation_analysis_metab <- run_mediation_analysis(
  data = combined_df_metab,
  input_cols = 5:6,
  output_cols = 7:8,
  mediator_cols = 9:264,
  covar_cols = 2:4,
  numOfCores = 36
)
write.csv(mediation_analysis_metab,"metab_mediation.csv")
mediation_analysis_kegg <- run_mediation_analysis(
  data = combined_df_kegg,
  input_cols = 5:6,
  output_cols = 7:8,
  mediator_cols = 9:2516,
  covar_cols = 2:4,
  numOfCores = 36
)
write.csv(mediation_analysis_kegg,"kegg_mediation.csv")
mediation_analysis_gmm <- run_mediation_analysis(
  data = combined_df_gmm,
  input_cols = 5:6,
  output_cols = 7:8,
  mediator_cols = 9:22,
  covar_cols = 2:4,
  numOfCores = 36
)
write.csv(mediation_analysis_gmm,"gmm_mediation.csv")
mediation_analysis_mgs <- run_mediation_analysis(
  data = combined_df_mgs,
  input_cols = 5:6,
  output_cols = 7:8,
  mediator_cols = 9:261,
  covar_cols = 2:4,
  numOfCores = 36
)
write.csv(mediation_analysis_mgs,"mgs_mediation.csv")
mediation_analysis <- rbind(mediation_analysis_metab,mediation_analysis_kegg,mediation_analysis_gmm,mediation_analysis_mgs)
# Flag top mediators
flag_top_mediators <- function(results_df, p_thresh = 0.05, top_prop = 0.1) {
  results_df %>%
    dplyr::filter(p_adj_boot < p_thresh) %>%
    dplyr::arrange(p_adj_boot, desc(abs(acme))) %>%
    dplyr::slice_head(n = ceiling(n() * top_prop))
}
top_hits <- flag_top_mediators(mediation_analysis,top_prop = 0.25)
# Write the results to a single CSV file
write.csv(mediation_analysis, "mediation_results.csv", row.names = FALSE)
