library(readr)
library(vegan)
library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(lefser)
library(mia)
library(ggplot2)
library(stringr)
library(LOCOM)

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
#ANCOM-BC2
perform_ancombc2_analysis <- function(otu_data, meta_data, group) {
  # Ensure metadata columns are correctly typed
  meta_data$AGE <- as.numeric(meta_data$AGE)
  meta_data$GENDER <- as.factor(meta_data$GENDER)
  meta_data$CENTER_C <- as.factor(meta_data$CENTER_C)
  
  if (!group %in% colnames(meta_data)) {
    stop(paste("Group variable", group, "not found in metadata."))
  }
  meta_data[[group]] <- as.factor(meta_data[[group]])
  
  # Remove samples with missing values
  covariates <- c("AGE", "GENDER", "CENTER_C", group)
  meta_data <- meta_data[complete.cases(meta_data[, covariates]), ]
  
  # Check for low sample size per group
  group_counts <- table(meta_data[[group]])
  if (any(group_counts < 5)) {
    warning("Some groups have fewer than 5 samples. Variance estimation may be unstable.")
  }
  
  # Drop CENTER_C if it has only one level
  if (length(unique(meta_data$CENTER_C)) < 2) {
    warning("CENTER_C has only one level and will be excluded from the model.")
    fix_formula <- paste("AGE + GENDER +", group)
  } else {
    fix_formula <- paste("AGE + GENDER + CENTER_C +", group)
  }
  
  # Match OTU data to metadata
  otu_data <- otu_data[match(meta_data$SampleID, rownames(otu_data)), ]
  
  # Create phyloseq object
  otu_mat <- as.matrix(otu_data)
  otu <- otu_table(otu_mat, taxa_are_rows = FALSE)
  meta <- sample_data(meta_data)
  rownames(meta) <- meta_data$SampleID
  ps_obj <- phyloseq(otu, meta)
  
  # Run ANCOM-BC2 with safe defaults
  output <- tryCatch({
    ancombc2(
      data = ps_obj,
      fix_formula = fix_formula,
      rand_formula = NULL,
      p_adj_method = "BH",
      pseudo_sens = TRUE,
      prv_cut = 0.01,
      lib_cut = 0,
      s0_perc = 0.05,
      group = group,
      struc_zero = FALSE,
      neg_lb = FALSE,
      alpha = 0.05,
      n_cl = 2,
      verbose = TRUE,
      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,
      iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = NULL,
      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
      trend_control = NULL
    )
  }, error = function(e) {
    warning("ANCOM-BC2 estimation failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(output)) {
    return(data.frame())
  }
  
  # Extract results
  res_prim <- output$res
  cat("Available covariates in beta:\n")
  print(names(res_prim$beta))
  
  cat("Available covariates in q_val:\n")
  print(names(res_prim$q_val))
  
  # Check if results are empty
  if (length(res_prim$beta) == 0 || length(res_prim$q_val) == 0) {
    warning("No taxa passed filtering or estimation failed. Returning empty result.")
    return(data.frame())
  }
  
  beta_group <- NULL
  qval_group <- NULL
  if (group %in% names(res_prim$beta)) {
    beta_group <- res_prim$beta[[group]]
    qval_group <- res_prim$q_val[[group]]
    cat("Preview of q-values for group:\n")
    print(head(qval_group))
  } else {
    warning(paste("Group", group, "not found in ANCOM-BC2 results."))
  }
  
  # Create results dataframe
  diff_abundant_df <- data.frame(
    Taxon = rownames(res_prim$beta),
    Beta_AGE = res_prim$beta$AGE,
    Beta_GENDER = res_prim$beta$GENDER,
    Beta_CENTER_C = if ("CENTER_C" %in% names(res_prim$beta)) res_prim$beta$CENTER_C else NA,
    Beta_GROUP = beta_group,
    Qval_AGE = res_prim$q_val$AGE,
    Qval_GENDER = res_prim$q_val$GENDER,
    Qval_CENTER_C = if ("CENTER_C" %in% names(res_prim$q_val)) res_prim$q_val$CENTER_C else NA,
    Qval_GROUP = qval_group
  )
  
  # Filter significant taxa
  significant_df <- diff_abundant_df[
    (diff_abundant_df$Qval_AGE < 0.05) &
      (diff_abundant_df$Qval_GENDER < 0.05) &
      (is.na(diff_abundant_df$Qval_CENTER_C) | diff_abundant_df$Qval_CENTER_C < 0.05) &
      (is.na(diff_abundant_df$Qval_GROUP) | diff_abundant_df$Qval_GROUP < 0.05),
  ]
  
  return(significant_df)
}
drug_non_drug <- perform_ancombc2_analysis(otu_data = drug_non_drug_mg,meta_data = drug_non_drug_md,group="Group")
drug_non_drug_0 <- perform_ancombc2_analysis(otu_data = drug_non_drug_0_mg,meta_data=drug_non_drug_0_md,group="Group")
run_lefse_from_counts <- function(otu_data, meta_data, group, reference_level = NULL) {
  library(SummarizedExperiment)
  library(lefser)
  
  # Ensure SampleID is rownames
  rownames(meta_data) <- meta_data$SampleID
  
  # Match samples
  common_samples <- intersect(rownames(meta_data), rownames(otu_data))
  otu_data <- otu_data[common_samples, , drop = FALSE]
  meta_data <- meta_data[common_samples, , drop = FALSE]
  
  # Relevel group if needed
  meta_data[[group]] <- as.factor(meta_data[[group]])
  if (!is.null(reference_level)) {
    meta_data[[group]] <- relevel(meta_data[[group]], ref = reference_level)
  }
  
  # Transpose OTU data: features as rows, samples as columns
  otu_data <- t(otu_data)
  colnames(otu_data) <- rownames(meta_data)
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(otu_data)),
    colData = S4Vectors::DataFrame(meta_data)
  )
  
  # Keep only terminal nodes
  terminal_features <- lefser::get_terminal_nodes(rownames(se))
  se <- se[terminal_features, ]
  
  # Compute relative abundances
  se_relab <- lefser::relativeAb(se)
  
  # Run LEfSe
  lefse_result <- lefser::lefser(
    relab = se_relab,
    classCol = group,
    kruskal.threshold = 0.05,
    wilcox.threshold = 0.05,
    lda.threshold = 2.0,
    assay = 1L,
    trim.names = FALSE,
    checkAbundances = TRUE,
    method= "BH",
  )
  return(lefse_result)
}
drug_non_drug_lefse <- run_lefse_from_counts(otu_data = drug_non_drug_mg,meta_data = drug_non_drug_md,group="Group")
drug_non_drug_0_lefse <- run_lefse_from_counts(otu_data = drug_non_drug_0_mg,meta_data=drug_non_drug_0_md,group="Group")
drug_non_drug_0_lefse <- dplyr::left_join(drug_non_drug_0_lefse,metagenome_annotation,by=c("features"="VariableID"))
drug_non_drug_lefse <- dplyr::left_join(drug_non_drug_lefse,metagenome_annotation,by=c("features"="VariableID"))

#Plot Function
plot_lda_scores <- function(df, df_name = "lda_plot") {
  df <- df %>%
    mutate(Species_Name = str_extract(Notes, "[^;]+$")) %>%
    filter(scores >= 2.5 | scores <= -2.5) %>%
    mutate(hjust_val = ifelse(scores > 0, 1, 0))
  
  p <- ggplot(df, aes(x = scores, y = reorder(DisplayName, scores))) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_text(aes(label = ifelse(!is.na(Species_Name), Species_Name, ""),
                  hjust = hjust_val),
              color = "black") +
    labs(title = "LDA Scores from LEfSe (|score| ≥ 2.5)",
         x = "LDA Score",
         y = "CAG ID") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank())
  
  print(p)
  
  # Save the plot
  ggsave(filename = paste0(df_name, ".png"), plot = p, width = 8, height = 10, dpi = 600)
}
drug_non_drug_lefse_plot <- plot_lda_scores(drug_non_drug_lefse,df_name = "drug_non_drug_lefse")
drug_non_drug_0_lefse_plot <- plot_lda_scores(drug_non_drug_0_lefse,df_name = "drug_non_drug_0_lefse")
drug_non_drug_diff_plot <- plot_lda_scores(drug_non_drug_diff,df_name = "drug_non_drug_diff")

#LOCOM
run_locom_analysis <- function(otu_data, meta_data) {
  # Ensure categorical variables are factors
  meta_data$GENDER <- as.factor(meta_data$GENDER)
  meta_data$CENTER_C <- as.factor(meta_data$CENTER_C)
  meta_data$Group <- as.factor(meta_data$Group)
  
  # Match OTU data to metadata
  otu_data <- otu_data[match(meta_data$SampleID, rownames(otu_data)), ]
  rownames(meta_data) <- meta_data$SampleID
  meta_data$SampleID <- NULL
  
  # Define Y (trait of interest) as binary (0/1)
  Y <- ifelse(meta_data$Group == levels(meta_data$Group)[1], 0, 1)
  
  # Create covariate matrix C safely
  C_formula <- ~ GENDER + AGE + CENTER_C
  C_model <- model.matrix(C_formula, data = meta_data)
  C <- C_model[, -1, drop = FALSE]  # Drop intercept
  storage.mode(C) <- "numeric"      # Ensure numeric type
  
  # Run LOCOM
  result <- locom(
    otu.table = otu_data,
    Y = Y,
    C = C,
    n.cores = 4,
    fdr.nominal = 0.05,
    filter.thresh = 0.05
  )
  
  # Format output
  taxa <- names(result$q.otu)
  output_df <- data.frame(
    Taxon = taxa,
    P_value = result$p.otu[taxa],
    Q_value = result$q.otu[taxa],
    Detected = taxa %in% result$detected.otu,
    Reference_Taxon = result$ref.tax
  )
  
  return(output_df)
}
drug_non_drug_locom <- run_locom_analysis(otu_data = drug_non_drug_mg,meta_data = drug_non_drug_md)
