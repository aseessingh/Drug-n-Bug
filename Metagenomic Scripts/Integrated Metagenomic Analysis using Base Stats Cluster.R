library(readr)
library(tidyverse)
library(stringr)

#Input Parameters
clustered_data <- as.data.frame(read_csv("meta_clustering.csv"))
metagenome <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv"))
drug_dose <- as.data.frame(read_csv("MC_Drug_doses.csv"))
metagenome_annotation <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.variables.v6.r.csv"))

#Preprocessing
metagenome$MGSampleID <- gsub("M0_","",metagenome$MGSampleID)
metagenome_annotation <- metagenome_annotation %>% dplyr::select(VariableID,DisplayName,Notes)
log_transform_normalize <- function(df) {
  # Exclude the first column
  data_to_transform <- df[, -1]
  
  # Apply log transformation (adding 1 to avoid log(0))
  log_transformed <- log(data_to_transform + 1)
  
  # Min-max normalization function
  normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Apply normalization
  normalized <- as.data.frame(lapply(log_transformed, normalize))
  
  # Combine the first column back with the transformed data
  result <- cbind(df[, 1, drop = FALSE], normalized)
  
  return(result)
}
metagenome <- log_transform_normalize(metagenome)
cluster_2 <- clustered_data %>% filter(Cluster == 2) %>% select(SampleID,Cluster,GENDER,CENTER_C,AGE)
cluster_1 <- clustered_data %>% filter (Cluster == 1) %>% select(SampleID,Cluster,GENDER,CENTER_C,AGE)
metagenome_1 <- metagenome %>% filter(MGSampleID %in% cluster_1$SampleID)
metagenome_2 <- metagenome %>% filter(MGSampleID %in% cluster_2$SampleID)
rownames_function <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}
metagenome_1 <- rownames_function(metagenome_1)
metagenome_2 <- rownames_function(metagenome_2)
metadata <- rbind(cluster_1,cluster_2)
metagenome_net <- rbind(metagenome_1,metagenome_2)
combined_analysis <- function(combined_df, metadata_df, annotation_df, alpha = 0.05,
                              covariates = c("AGE", "GENDER", "CENTER_C"), name_key = NULL) {
  library(dplyr)
  library(stats)
  library(broom)
  
  combined_df$SampleID <- rownames(combined_df)
  if (!"SampleID" %in% colnames(metadata_df)) {
    metadata_df$SampleID <- rownames(metadata_df)
  }
  
  merged_df <- merge(metadata_df, combined_df, by = "SampleID")
  if (!"Cluster" %in% colnames(merged_df)) stop("Cluster column not found in metadata.")
  merged_df$Cluster <- factor(merged_df$Cluster)
  
  exclude_cols <- c("SampleID", "Cluster", covariates)
  feature_cols <- setdiff(colnames(merged_df), exclude_cols)
  
  # Step 1: Mann-Whitney U Test
  mw_results <- lapply(feature_cols, function(feature) {
    values <- merged_df[[feature]]
    cluster <- merged_df$Cluster
    if (length(unique(cluster)) != 2) return(NULL)
    
    cluster1_values <- values[cluster == levels(cluster)[1]]
    cluster2_values <- values[cluster == levels(cluster)[2]]
    cluster1_values <- cluster1_values[!is.na(cluster1_values)]
    cluster2_values <- cluster2_values[!is.na(cluster2_values)]
    if (length(cluster1_values) < 2 || length(cluster2_values) < 2) return(NULL)
    
    test <- tryCatch(wilcox.test(cluster1_values, cluster2_values, exact = FALSE), error = function(e) NULL)
    if (is.null(test)) return(NULL)
    
    median1 <- median(cluster1_values)
    median2 <- median(cluster2_values)
    higher_cluster <- ifelse(median1 > median2, levels(cluster)[1], levels(cluster)[2])
    n1 <- length(cluster1_values)
    n2 <- length(cluster2_values)
    u <- test$statistic
    rbc <- (2 * u) / (n1 * n2) - 1
    
    data.frame(
      Feature = feature,
      MW_U_Statistic = as.numeric(u),
      MW_P_Value = test$p.value,
      MW_Effect_Size = rbc,
      MW_Higher_Cluster = higher_cluster
    )
  })
  
  mw_df <- do.call(rbind, mw_results)
  if (is.null(mw_df) || nrow(mw_df) == 0) return(NULL)
  
  mw_df$MW_Adjusted_P <- p.adjust(mw_df$MW_P_Value, method = "BH")
  mw_df <- mw_df[mw_df$MW_Adjusted_P < alpha, ]
  if (nrow(mw_df) == 0) return(NULL)
  
  relevant_features <- mw_df$Feature
  
  # Step 2: GLM
  glm_results <- lapply(relevant_features, function(feature) {
    formula <- as.formula(paste0("`", feature, "` ~ Cluster + ", paste(covariates, collapse = " + ")))
    model <- tryCatch(glm(formula, data = merged_df, family = gaussian()), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    coef_name <- grep("^Cluster", rownames(summary(model)$coefficients), value = TRUE)[1]
    if (!is.na(coef_name)) {
      coef <- summary(model)$coefficients[coef_name, "Estimate"]
      p_value <- summary(model)$coefficients[coef_name, "Pr(>|t|)"]
    } else {
      coef <- NA
      p_value <- NA
    }
    data.frame(Feature = feature, GLM_Beta = coef, GLM_P_Value = p_value)
  })
  
  glm_df <- do.call(rbind, glm_results)
  glm_df$GLM_Adjusted_P <- p.adjust(glm_df$GLM_P_Value, method = "BH")
  
  # Step 3: Logistic Regression
  logit_results <- lapply(relevant_features, function(feature) {
    formula <- as.formula(paste0("Cluster ~ `", feature, "` + ", paste(covariates, collapse = " + ")))
    model <- tryCatch(glm(formula, data = merged_df, family = binomial()), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    coef_name <- grep(paste0("^`?", feature, "`?$"), rownames(summary(model)$coefficients), value = TRUE)[1]
    if (!is.na(coef_name)) {
      coef <- summary(model)$coefficients[coef_name, "Estimate"]
      p_value <- summary(model)$coefficients[coef_name, "Pr(>|z|)"]
      odds_ratio <- exp(coef)
    } else {
      odds_ratio <- NA
      p_value <- NA
    }
    data.frame(Feature = feature, Logit_OR = odds_ratio, Logit_P_Value = p_value)
  })
  
  logit_df <- do.call(rbind, logit_results)
  logit_df$Logit_Adjusted_P <- p.adjust(logit_df$Logit_P_Value, method = "BH")
  
  # Step 4: Merge and Annotate
  result_df <- mw_df %>%
    inner_join(glm_df, by = "Feature") %>%
    inner_join(logit_df, by = "Feature") %>%
    mutate(
      MW_Significant = MW_Adjusted_P < alpha,
      GLM_Significant = GLM_Adjusted_P < alpha,
      Logit_Significant = Logit_Adjusted_P < alpha,
      Significance_Score = rowSums(across(c(MW_Significant, GLM_Significant, Logit_Significant))),
      Rank = rank(-Significance_Score, ties.method = "min")
    ) %>%
    arrange(desc(Significance_Score), MW_Adjusted_P)
  
  if (!is.null(name_key)) {
    result_df <- result_df %>%
      left_join(name_key, by = c("Feature" = "cleaned")) %>%
      mutate(Feature = ifelse(!is.na(original), original, Feature)) %>%
      select(-original)
  }
  
  result_df <- result_df %>%
    left_join(annotation_df, by = c("Feature" = "VariableID"))
  
  return(result_df)
}
results <- combined_analysis(combined_df = metagenome,alpha = 0.05,metadata_df = metadata,annotation_df=metagenome_annotation)
write.csv(results,"cluster_data_metagenome.csv")
