library(readr)
library(tidyverse)
library(metadeconfoundR)
library(stringr)
library(ggplot2)
library(ggrepel)

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

#MetadeconfoundR
drug_non_drug_md <- drug_non_drug_md %>% distinct(SampleID,.keep_all = TRUE)
rownames(drug_non_drug_md) <- drug_non_drug_md[,1]
drug_non_drug_md[,1] <- NULL
drug_non_drug_md <- drug_non_drug_md[match(rownames(drug_non_drug_mg), rownames(drug_non_drug_md)), ]
drug_non_drug_metadeconfound_run<- MetaDeconfound(
  featureMat = drug_non_drug_mg,
  metaMat = drug_non_drug_md,
  typeCategorical = c("GENDER","CENTER_C"),
  typeContinuous = c("AGE"),
  QCutoff = 0.05,
  DCutoff = 0.2,
  nnodes = 16,
  adjustMethod = "fdr"
)
drug_non_drug_metadeconfound <- as.data.frame(cbind(drug_non_drug_metadeconfound_run$Ps,drug_non_drug_metadeconfound_run$Qs,drug_non_drug_metadeconfound_run$Ds,drug_non_drug_metadeconfound_run$status))
colnames(drug_non_drug_metadeconfound) <- make.unique(colnames(drug_non_drug_metadeconfound))
drug_non_drug_metadeconfound <- drug_non_drug_metadeconfound %>%
  filter(Group.3 %in% c("OK_nc", "OK_sd")) %>% select (Group,Group.1,Group.2) %>% rename(p = Group, q = Group.1 , Rho = Group.2)
drug_non_drug_metadeconfound$Variable <- rownames(drug_non_drug_metadeconfound)
drug_non_drug_metadeconfound <- left_join(drug_non_drug_metadeconfound,metagenome_annotation,by=c("Variable"="VariableID"))
rownames(drug_non_drug_metadeconfound) <- NULL
write.csv(drug_non_drug_metadeconfound,"drug_non_drug_metadeconfound.csv")
#
drug_non_drug_0_md <- drug_non_drug_0_md %>% distinct(SampleID,.keep_all = TRUE)
rownames(drug_non_drug_0_md) <- drug_non_drug_0_md[,1]
drug_non_drug_0_md[,1] <- NULL
drug_non_drug_0_md <- drug_non_drug_0_md[match(rownames(drug_non_drug_0_mg), rownames(drug_non_drug_0_md)), ]
drug_non_drug_0_metadeconfound_run<- MetaDeconfound(
  featureMat = drug_non_drug_0_mg,
  metaMat = drug_non_drug_0_md,
  typeCategorical = c("GENDER","CENTER_C"),
  typeContinuous = c("AGE"),
  QCutoff = 0.05,
  DCutoff = 0.2,
  nnodes = 16,
  adjustMethod = "fdr"
)
drug_non_drug_0_metadeconfound <- as.data.frame(cbind(drug_non_drug_0_metadeconfound_run$Ps,drug_non_drug_0_metadeconfound_run$Qs,drug_non_drug_0_metadeconfound_run$Ds,drug_non_drug_0_metadeconfound_run$status))
colnames(drug_non_drug_0_metadeconfound) <- make.unique(colnames(drug_non_drug_0_metadeconfound))
drug_non_drug_0_metadeconfound <- drug_non_drug_0_metadeconfound %>%
  filter(Group.3 %in% c("OK_nc", "OK_sd")) %>% select (Group,Group.1,Group.2) %>% rename(p = Group, q = Group.1 , Rho = Group.2)
drug_non_drug_0_metadeconfound$Variable <- rownames(drug_non_drug_0_metadeconfound)
drug_non_drug_0_metadeconfound <- left_join(drug_non_drug_0_metadeconfound,metagenome_annotation,by=c("Variable"="VariableID"))
rownames(drug_non_drug_0_metadeconfound) <- NULL
write.csv(drug_non_drug_0_metadeconfound,"drug_non_drug_0_metadeconfound.csv")

#MWU-GLM-LR
combined_analysis <- function(combined_df, metadata_df, annotation_df, alpha = 0.05, covariates = c("AGE", "GENDER", "CENTER_C")) {
  library(dplyr)
  library(stats)
  library(broom)
  
  # Add SampleID as a column
  combined_df$SampleID <- rownames(combined_df)
  metadata_df$SampleID <- as.character(metadata_df$SampleID)
  
  # Merge metadata and features
  merged_df <- merge(metadata_df, combined_df, by = "SampleID")
  
  # Ensure Group is a factor
  if (!"Group" %in% colnames(merged_df)) stop("Group column not found in metadata.")
  merged_df$Group <- factor(merged_df$Group)
  
  # Identify feature columns (exclude SampleID and covariates)
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
    left_join(annotation_df, by = c("Feature" = "VariableID"))
  
  return(result_df)
}
drug_non_drug_net <- combined_analysis(combined_df = drug_non_drug_mg,alpha = 0.05,metadata_df = drug_non_drug_md,annotation_df=metagenome_annotation)
drug_non_drug_0_net <- combined_analysis(combined_df = drug_non_drug_0_mg,alpha = 0.05,metadata_df = drug_non_drug_0_md,annotation_df=metagenome_annotation)
write.csv(drug_non_drug_net,"drug_non_drug_base_stats.csv")
write.csv(drug_non_drug_0_net,"drug_non_drug_0_base_stats.csv")

#Plots
# Prepare data
results <- drug_non_drug_net
results <- results %>% select(GLM_Beta, Logit_OR, DisplayName, Notes)
results$Logit_OR <- log10(results$Logit_OR)
results$GLM_Beta <- 10^(results$GLM_Beta)
results_metad <- drug_non_drug_metadeconfound

# Scatter plot function
plot_glm_vs_logit <- function(data,
                              x_col = "GLM_Beta",
                              y_col = "Logit_OR",
                              display_col = "DisplayName",
                              notes_col = "Notes") {
  data <- data %>%
    mutate(
      Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
      Species = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 7],
      Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum),
      Species = ifelse(is.na(Species) | Species == "", !!sym(display_col), Species),
      Label = paste0(!!sym(display_col), "\n", Species)
    )
  
  x_thresh <- quantile(abs(data[[x_col]]), 0.95, na.rm = TRUE)
  y_thresh <- quantile(abs(data[[y_col]]), 0.95, na.rm = TRUE)
  
  data_top <- data %>%
    filter(abs(!!sym(x_col)) >= x_thresh | abs(!!sym(y_col)) >= y_thresh)
  
  ggplot(data_top, aes_string(x = x_col, y = y_col, label = "Label", color = "Phylum")) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, color = "gray40", linetype = "dotted") +
    geom_text_repel(size = 3, max.overlaps = Inf) +
    labs(
      x = "10^GLM Beta Coefficient",
      y = "log10(Odds Ratio)",
      title = "Scatter Plot of Beta Coefficient from GLM and Logistic Regression",
      color = "Phylum"
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}

# Bar plot function
plot_rho_bar <- function(data,
                         rho_col = "Rho",
                         display_col = "DisplayName",
                         notes_col = "Notes") {
  data <- data %>%
    mutate(
      Phylum = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 2],
      Species = str_split(!!sym(notes_col), ";", simplify = TRUE)[, 7],
      Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum),
      Species = ifelse(is.na(Species) | Species == "", !!sym(display_col), Species)
    )
  
  rho_thresh <- quantile(abs(data[[rho_col]]), 0.7, na.rm = TRUE)
  
  data_top <- data %>%
    filter(abs(.data[[rho_col]]) >= rho_thresh) %>%
    arrange(.data[[rho_col]]) %>%
    mutate(!!display_col := factor(.data[[display_col]], levels = .data[[display_col]]))
  
  ggplot(data_top, aes(x = .data[[rho_col]], y = .data[[display_col]], fill = Phylum)) +
    geom_col() +
    geom_text(aes(label = Species, hjust = ifelse(.data[[rho_col]] > 0, -0.1, 1.1)), size = 3) +
    labs(
      x = "Rho",
      y = "CAGID",
      title = "Correlation Coefficient Plot (Top 30%)",
      fill = "Phylum"
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}

# Generate plots
plot1 <- plot_glm_vs_logit(data = results)
plot2 <- plot_rho_bar(data = results_metad)

# Save plots
ggsave(filename = "base_stats.png", plot = plot1, dpi = 600, width = 14, height = 10)
ggsave(filename = "metadeconfoundR.png", plot = plot2, dpi = 600, width = 14, height = 10)
