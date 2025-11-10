library(readr)
library(tidyverse)
library(metadeconfoundR)

#Input Parameters
drug_intake <- as.data.frame(read_csv("metadata3.csv"))
non_drug_intake <- as.data.frame(read_csv("metadata4.csv"))
metagenome <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv"))
drug_dose <- as.data.frame(read_csv("MC_Drug_doses.csv"))
drug <- c("metformin")
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

#KW-MWU-GLM-LogR
combined_analysis <- function(group1_df, group2_df, metadata_df, annotation_df, alpha = 0.05, covariates = c("AGE", "GENDER", "CENTER_C")) {
  library(dplyr)
  library(stats)
  library(broom)
  
  # Step 1: Kruskal-Wallis Test
  combined_df <- rbind(group1_df, group2_df)
  kw_results <- lapply(names(combined_df), function(feature) {
    group1_values <- group1_df[[feature]]
    group2_values <- group2_df[[feature]]
    test <- tryCatch({
      kruskal.test(list(group1_values, group2_values))
    }, error = function(e) NULL)
    
    if (!is.null(test)) {
      data.frame(Feature = feature, KW_Statistic = test$statistic, KW_P_Value = test$p.value)
    } else {
      data.frame(Feature = feature, KW_Statistic = NA, KW_P_Value = NA)
    }
  })
  kw_df <- do.call(rbind, kw_results)
  kw_df$KW_Adjusted_P <- p.adjust(kw_df$KW_P_Value, method = "BH")
  kw_df <- kw_df[!is.na(kw_df$KW_Adjusted_P) & kw_df$KW_Adjusted_P < alpha, ]
  
  # Step 2: Mann-Whitney U Test
  mw_results <- lapply(kw_df$Feature, function(feature) {
    group1_values <- group1_df[[feature]]
    group2_values <- group2_df[[feature]]
    test <- wilcox.test(group1_values, group2_values, exact = FALSE)
    median1 <- median(group1_values, na.rm = TRUE)
    median2 <- median(group2_values, na.rm = TRUE)
    higher_group <- ifelse(median1 > median2, "Group 1", "Group 2")
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
  mw_df$MW_Adjusted_P <- p.adjust(mw_df$MW_P_Value, method = "BH")
  mw_df <- mw_df[mw_df$MW_Adjusted_P < alpha, ]
  
  # Step 3: Prepare merged data
  relevant_features <- mw_df$Feature
  metagenome_df <- rbind(group1_df, group2_df)[, relevant_features, drop = FALSE]
  metagenome_df$SampleID <- rownames(metagenome_df)
  metadata_df$SampleID <- rownames(metadata_df)
  merged_df <- merge(metadata_df, metagenome_df, by = "SampleID")
  merged_df$Group <- factor(merged_df$Group)
  
  # Step 4: GLM
  glm_results <- lapply(relevant_features, function(feature) {
    formula <- as.formula(paste0("`", feature, "` ~ Group + ", paste(covariates, collapse = " + ")))
    model <- glm(formula, data = merged_df, family = gaussian())
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
  
  # Step 5: Logistic Regression
  logit_results <- lapply(relevant_features, function(feature) {
    formula <- as.formula(paste0("Group ~ `", feature, "` + ", paste(covariates, collapse = " + ")))
    model <- glm(formula, data = merged_df, family = binomial())
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
  
  # Step 6: Merge all results
  result_df <- kw_df %>%
    inner_join(mw_df, by = "Feature") %>%
    inner_join(glm_df, by = "Feature") %>%
    inner_join(logit_df, by = "Feature")
  
  # Step 7: Annotate
  final_df <- result_df %>%
    left_join(annotation_df, by = c("Feature" = "VariableID"))
  
  return(final_df)
}
drug_non_drug_net <- combined_analysis(group1_df = metagenome_drug,group2_df = metagenome_non_drug,alpha = 0.05,metadata_df = drug_non_drug_md,annotation_df=metagenome_annotation)
drug_non_drug_0_net <- combined_analysis(group1_df = metagenome_drug,group2_df = metagenome_non_drug_0,alpha = 0.05,metadata_df = drug_non_drug_0_md,annotation_df=metagenome_annotation)
write.csv(drug_non_drug_net,"drug_non_drug_base_stats.csv")
write.csv(drug_non_drug_0_net,"drug_non_drug_0_base_stats.csv")
