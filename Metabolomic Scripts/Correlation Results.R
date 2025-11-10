library(readr)
library(dplyr)
library(ggplot2)

clustered_data <- as.data.frame(read_csv("2 var clustering.csv",show_col_types = FALSE))
metadata <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
clustered_data <- clustered_data %>% select(SampleID,Cluster)
drug_dose <- as.data.frame(read_csv("MC_Drug_doses.csv",show_col_types = FALSE))
check_cluster_medication_correlation <- function(metadata, clustered_data, drug_dose = NULL) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stats)
  library(patchwork)
  
  # Lookup table
  lookup <- data.frame(
    Variable = c("METFORMIN_C", "SU_C", "DPPIV_C", "INSULIN_C", "ACARBOSE_C", "BOLUS_C", "GLITAZONE_C", "GLP_1_C", "SGLT2_C",
                 "AT2_INHIB_C", "THIAZIDIQUE_C", "CTRBLOCK_C", "ACE_INHIB_C", "BETA_BL_C", "CA2_CBL_C", "DIURANSE_C", "ANTAGRENINE_C", "NITRATE_C", "KSPDIUR_C",
                 "STATINE_C", "EZETIMIBE_C", "FIBRATE_C", "REDRICEOMEGA_C",
                 "AMIODARON_C", "DIGITALIS_C", "OTHER_ANTIARYTHMIC_C",
                 "ASA_C", "CLOPIDOGREL_C", "CUMARINE_C", "HEPARINE_C", "ANTITHOMBO_C",
                 "PPI_AND_PPI_RELATED_C", "PPI_C", "PPI_RELATED_C",
                 "XANTOXINH_C"),
    Label = c("Metformin", "Sulfonylureas", "DPP-4 Inhibitors", "Insulin", "Acarbose", "Bolus Insulin", "PPAR-gamma Inhibitors", "GLP-1 Agonists", "SGLT2 Inhibitors",
              "AT2 Inhibitors", "Thiazides", "Calcium Channel Blockers", "ACE Inhibitors", "Beta Blockers", "Calcium Channel Blockers", "Diuretics", "Angiotensin II Receptor Blockers", "Nitrates", "Potassium-sparing Diuretics",
              "Statins", "Ezetimibe", "Fibrates", "Red Yeast Rice & Omega-3",
              "Amiodarone", "Digitalis", "Other Antiarrhythmics",
              "Aspirin", "Clopidogrel", "Coumarins", "Heparin", "Antithrombotics",
              "PPI and Related", "Proton Pump Inhibitors", "PPI Related",
              "Xanthine Oxidase Inhibitors"),
    stringsAsFactors = FALSE
  )
  
  categories <- list(
    antidiabetic = lookup$Variable[1:9],
    antihypert = lookup$Variable[10:19],
    antilipid = lookup$Variable[20:23],
    antiarhytmetic = lookup$Variable[24:26],
    anticoagulant = lookup$Variable[27:31],
    PPI = lookup$Variable[32:34],
    antigout = lookup$Variable[35]
  )
  
  metadata_joined <- metadata %>%
    inner_join(clustered_data, by = "SampleID")
  
  results <- data.frame(Category = character(), Chi2 = numeric(), P_value = numeric(),
                        Cluster = character(), Proportion = numeric(), Group = character(), stringsAsFactors = FALSE)
  
  for (cat in names(categories)) {
    meds <- categories[[cat]]
    meds_present <- meds[meds %in% colnames(metadata_joined)]
    
    if (length(meds_present) == 0) next
    
    if (length(meds_present) == 1) {
      metadata_joined[[cat]] <- metadata_joined[[meds_present]] > 0
    } else {
      metadata_joined[[cat]] <- rowSums(metadata_joined[, meds_present], na.rm = TRUE) > 0
    }
    
    tbl_cat <- table(metadata_joined$Cluster, metadata_joined[[cat]])
    if (all(dim(tbl_cat) > 1)) {
      test_cat <- suppressWarnings(chisq.test(tbl_cat))
      for (cl in rownames(tbl_cat)) {
        prop <- tbl_cat[cl, "TRUE"] / sum(tbl_cat[cl, ])
        results <- rbind(results, data.frame(Category = cat, Chi2 = test_cat$statistic,
                                             P_value = test_cat$p.value, Cluster = cl,
                                             Proportion = prop, Group = "Category"))
      }
    }
    
    for (med in meds_present) {
      tbl_med <- table(metadata_joined$Cluster, metadata_joined[[med]])
      if (all(dim(tbl_med) > 1)) {
        test_med <- suppressWarnings(chisq.test(tbl_med))
        for (cl in rownames(tbl_med)) {
          prop <- tbl_med[cl, "1"] / sum(tbl_med[cl, ])
          results <- rbind(results, data.frame(Category = med, Chi2 = test_med$statistic,
                                               P_value = test_med$p.value, Cluster = cl,
                                               Proportion = prop, Group = cat))
        }
      }
    }
  }
  
  results$FDR <- p.adjust(results$P_value, method = "fdr")
  results$Significance <- cut(results$FDR,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", ""),
                              right = TRUE)
  
  results <- results %>%
    group_by(Category) %>%
    filter(any(Proportion > 0)) %>%
    ungroup() %>%
    left_join(lookup, by = c("Category" = "Variable")) %>%
    mutate(Label = ifelse(is.na(Label), Category, Label))
  
  heatmap_plot <- ggplot(results, aes(x = Cluster, y = Label, fill = Proportion)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Significance), color = "black", size = 3, na.rm = TRUE) +
    scale_fill_gradient(low = "blue", high = "green", name = "Proportion") +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Cluster vs Medication Association", x = "Cluster", y = "Drug Class") +
    theme_minimal()
  
  correlation_plot <- NULL
  if (!is.null(drug_dose)) {
    dose_joined <- drug_dose %>%
      rename(SampleID = MCid) %>%
      inner_join(clustered_data, by = "SampleID")
    
    dose_cols <- grep("^DOSAGE_", colnames(dose_joined), value = TRUE)
    dose_data <- dose_joined[, c("Cluster", dose_cols)]
    dose_data$Cluster <- as.numeric(as.factor(dose_data$Cluster))
    
    valid_dose_cols <- dose_cols[sapply(dose_data[, dose_cols], function(x) sd(x, na.rm = TRUE) > 0)]
    
    cor_vals <- sapply(valid_dose_cols, function(col) suppressWarnings(
      cor(dose_data[[col]], dose_data$Cluster, method = "spearman", use = "pairwise.complete.obs")
    ))
    
    p_vals <- sapply(valid_dose_cols, function(col) suppressWarnings(
      cor.test(dose_data[[col]], dose_data$Cluster, method = "spearman")$p.value
    ))
    
    correlation_matrix <- data.frame(
      Drug = gsub("^DOSAGE_", "", valid_dose_cols),
      Spearman_Correlation = cor_vals,
      P_value = p_vals,
      FDR = p.adjust(p_vals, method = "fdr")
    ) %>%
      mutate(Variable = Drug) %>%
      left_join(lookup, by = "Variable") %>%
      mutate(Label = ifelse(is.na(Label), Drug, Label),
             Significance = cut(FDR,
                                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                labels = c("***", "**", "*", ""),
                                right = TRUE))
    
    correlation_plot <- ggplot(correlation_matrix, aes(x = "Cluster", y = Label, fill = Spearman_Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Significance), color = "black", size = 3, na.rm = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Spearman ρ") +
      labs(title = "Spearman Correlation of Drug Dosage with Cluster", x = "", y = "Drug Dosage") +
      theme_minimal()
  }
  
  combined_plot <- if (!is.null(correlation_plot)) {
    heatmap_plot + correlation_plot +
      plot_layout(ncol = 2, widths = c(1, 1)) +
      plot_annotation(title = "Cluster Medication and Dosage Correlation Overview")
  } else {
    heatmap_plot
  }
  
  return(list(
    results = results,
    correlation_matrix = if (exists("correlation_matrix")) correlation_matrix else NULL,
    plot = combined_plot
  ))
}
results <- check_cluster_medication_correlation(metadata=metadata,clustered_data=clustered_data,drug_dose = drug_dose)
ggsave("correl and chi.png",plot = results$plot,dpi=600,width=18,height=10)
write.csv(results$correlation_matrix,"correl_cluster.csv")
write.csv(results$results,"chi_cluster.csv")



