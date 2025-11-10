library(readr)
library(tidyverse)
library(stringr)

#Input Parameters
drug_intake <- as.data.frame(read_csv("metadata3.csv"))
non_drug_intake <- as.data.frame(read_csv("metadata4.csv"))
metagenome <- as.data.frame(read_csv("hub.function_adjusted.GMM.down.10000000.fpkm.mgsamples.v5.r.csv"))
drug_dose <- as.data.frame(read_csv("MC_Drug_doses.csv"))
drug <- c("metformin")
metagenome_annotation <- as.data.frame(read_csv("hub.function_adjusted.GMM.down.10000000.fpkm.variables.v5.r.csv"))
#Preprocessing
metagenome$MGSampleID <- gsub("M0_","",metagenome$MGSampleID)
metagenome_annotation <- metagenome_annotation %>% dplyr::select(VariableID,DisplayName,Name,HL1,HL2)
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
non_drug_intake <- non_drug_intake %>% filter(Status == 0)
drug_intake <- drug_intake %>% select(SampleID,AGE,GENDER,CENTER_C)
non_drug_intake <- dplyr::left_join(non_drug_intake,drug_dose,by=c("SampleID"="MCid"))
filter_non_drug_intake <- function(df, drug) {
  if (drug == "metformin") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_DPPIV_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_SU_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else if (drug == "su") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_DPPIV_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_METFORMIN_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else if (drug == "dppiv") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_SU_C == 0 & DOSAGE_METFORMIN_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else {
    stop("Invalid drug name")
  }
  
  return(non_drug_intake_0)
}
non_drug_intake_0 <- filter_non_drug_intake(df=non_drug_intake,drug=drug)
non_drug_intake <- non_drug_intake %>% select(SampleID,AGE,GENDER,CENTER_C)
metagenome_drug <- metagenome %>% filter(MGSampleID %in% drug_intake$SampleID)
metagenome_non_drug <- metagenome %>% filter(MGSampleID %in% non_drug_intake$SampleID)
metagenome_non_drug_0 <- metagenome %>% filter (MGSampleID %in% non_drug_intake_0$SampleID)
rownames_function <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}
metagenome_drug <- rownames_function(metagenome_drug)
metagenome_non_drug <- rownames_function(metagenome_non_drug)
metagenome_non_drug_0 <- rownames_function(metagenome_non_drug_0)
drug_intake$Group <- 1
non_drug_intake$Group <- 0
non_drug_intake_0$Group <- 0
drug_intake <- drug_intake %>% select(SampleID,Group, everything())
non_drug_intake <- non_drug_intake %>% select(SampleID,Group, everything())
non_drug_intake_0 <- non_drug_intake_0 %>% select(SampleID,Group, everything())
drug_non_drug_md <- rbind(drug_intake,non_drug_intake)
drug_non_drug_0_md <- rbind(drug_intake,non_drug_intake_0)
drug_non_drug_mg <- rbind(metagenome_drug,metagenome_non_drug)
drug_non_drug_0_mg <- rbind (metagenome_drug,metagenome_non_drug_0)

#MWU-GLM-LR
combined_analysis <- function(combined_df, metadata_df, annotation_df, alpha = 0.05,
                              covariates = c("AGE", "GENDER", "CENTER_C")) {
  library(dplyr)
  library(stats)
  library(broom)
  
  combined_df$SampleID <- rownames(combined_df)
  if (!"SampleID" %in% colnames(metadata_df)) {
    metadata_df$SampleID <- rownames(metadata_df)
  }
  
  merged_df <- merge(metadata_df, combined_df, by = "SampleID")
  if (!"Group" %in% colnames(merged_df)) stop("Group column not found in metadata.")
  merged_df$Group <- factor(merged_df$Group)
  
  exclude_cols <- c("SampleID", "Group", covariates)
  feature_cols <- setdiff(colnames(merged_df), exclude_cols)
  
  # Step 1: Mann-Whitney U Test
  mw_results <- lapply(feature_cols, function(feature) {
    values <- merged_df[[feature]]
    group <- merged_df$Group
    if (length(unique(group)) != 2) return(NULL)
    
    group1_values <- values[group == levels(group)[1]]
    group2_values <- values[group == levels(group)[2]]
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    if (length(group1_values) < 2 || length(group2_values) < 2) return(NULL)
    
    test <- tryCatch(wilcox.test(group1_values, group2_values, exact = FALSE), error = function(e) NULL)
    if (is.null(test)) return(NULL)
    
    median1 <- median(group1_values)
    median2 <- median(group2_values)
    higher_group <- ifelse(median1 > median2, levels(group)[1], levels(group)[2])
    n1 <- length(group1_values)
    n2 <- length(group2_values)
    u <- test$statistic
    rbc <- (2 * u) / (n1 * n2) - 1
    
    data.frame(
      Feature = feature,
      MW_U_Statistic = as.numeric(u),
      MW_P_Value = test$p.value,
      MW_Effect_Size = rbc,
      MW_Higher_Group = higher_group
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
    formula <- as.formula(paste0("`", feature, "` ~ Group + ", paste(covariates, collapse = " + ")))
    model <- tryCatch(glm(formula, data = merged_df, family = gaussian()), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    
    coef_name <- grep("^Group", rownames(summary(model)$coefficients), value = TRUE)[1]
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
    formula <- as.formula(paste0("Group ~ `", feature, "` + ", paste(covariates, collapse = " + ")))
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
    arrange(desc(Significance_Score), MW_Adjusted_P) %>%
    left_join(annotation_df, by = c("Feature" = "VariableID"))
  
  return(result_df)
}

drug_non_drug_net <- combined_analysis(combined_df = drug_non_drug_mg,alpha = 0.05,metadata_df = drug_non_drug_md,annotation_df=metagenome_annotation)
drug_non_drug_0_net <- combined_analysis(combined_df = drug_non_drug_0_mg,alpha = 0.05,metadata_df = drug_non_drug_0_md,annotation_df=metagenome_annotation)
write.csv(drug_non_drug_net,"drug_non_drug_base_mgs_stats.csv")
write.csv(drug_non_drug_0_net,"drug_non_drug_0_base_mgs_stats.csv")
