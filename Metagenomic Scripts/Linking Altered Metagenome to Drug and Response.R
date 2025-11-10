library(readr)
library(dplyr)

#Input
metadata <- as.data.frame(read_csv("metadata3.csv",show_col_types = FALSE))
group_3 <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
metabolite_for_AG <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))
metagenome <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv",show_col_types = FALSE))
results <- as.data.frame(read_csv("metformin vs non-metformin MGS.csv",show_col_types = FALSE))
dose_data <- c("DOSAGE_METFORMIN_C")
drug_metabolite <- c("metformin")

#Preprocessing
group_3 <- group_3 %>% select(SampleID,GLYCATHB)
metadata <- left_join(metadata,group_3,by=c("SampleID"))
metabolite_for_AG <- metabolite_for_AG %>% select(MGS_PID,`1,5-anhydroglucitol (1,5-AG)`) %>% rename (AG=`1,5-anhydroglucitol (1,5-AG)`)
metadata <- left_join(metadata,metabolite_for_AG,by=c("SampleID"="MGS_PID"))
metagenome <- metagenome %>% select(MGSampleID,results$Feature)
metagenome$MGSampleID <- gsub("M0_","",metagenome$MGSampleID)
metadata <- left_join(metadata,metagenome,by=c("SampleID"="MGSampleID"))

#Lm Function
lm_analysis <- function(df, metagenome_cols, drug_dose_col, drug_metabolite_col) {
  library(dplyr)
  library(broom)
  library(tibble)
  
  # Function to adjust p-values using FDR
  adjust_pvalues <- function(pvalues) {
    p.adjust(pvalues, method = "fdr")
  }
  
  perform_linear_regression <- function(combined_data, dependent_var) {
    all_p_values <- data.frame()
    
    # Ensure column names are unique
    colnames(combined_data) <- make.unique(colnames(combined_data))
    
    for (i in metagenome_cols) {  # Metagenomic data columns
      metabolite <- colnames(combined_data)[i]
      metabolite <- trimws(metabolite)  # Trim any extra spaces
      # Construct the formula correctly with quotes around the metabolite name
      formula <- as.formula(paste("`", metabolite, "` ~ `", dependent_var, "` + AGE + GENDER + CENTER_C", sep = ""))
      model <- lm(formula, data = combined_data)
      summary_model <- summary(model)
      
      p_value <- summary_model$coefficients[2, "Pr(>|t|)"]  # The dependent variable is always the second coefficient
      beta_coefficient <- summary_model$coefficients[2, "Estimate"]
      all_p_values <- rbind(all_p_values, data.frame(metabolite = metabolite, p_value = p_value, beta_coefficient = beta_coefficient))
    }
    
    # Apply FDR correction
    all_p_values$adjusted_p_value <- adjust_pvalues(all_p_values$p_value)
    
    # Extract significant results
    significant_results <- all_p_values %>% filter(adjusted_p_value < 0.05)
    
    return(significant_results)
  }
  
  # Perform linear regression for dose and metabolite data
  significant_results_dose <- perform_linear_regression(df, drug_dose_col)
  significant_results_metabolite <- perform_linear_regression(df, drug_metabolite_col)
  
  # Perform linear regression for response1 (GLYCATHB) and response2 (AG)
  significant_results_response1 <- perform_linear_regression(df, "GLYCATHB")
  significant_results_response2 <- perform_linear_regression(df, "AG")
  
  # Export results to CSV files
  write.csv(significant_results_dose, paste0(drug_dose_col, ".csv"), row.names = FALSE)
  write.csv(significant_results_metabolite, paste0(drug_metabolite_col, ".csv"), row.names = FALSE)
  write.csv(significant_results_response1, "GLYCATHB.csv", row.names = FALSE)
  write.csv(significant_results_response2, "AG.csv", row.names = FALSE)
  
  return(list(
    dose_data_results = significant_results_dose,
    metabolite_data_results = significant_results_metabolite,
    response1_results = significant_results_response1,
    response2_results = significant_results_response2
  ))
}
lm <- lm_analysis(
  df = metadata,
  metagenome_cols = 8:17,
  drug_dose_col = dose_data,
  drug_metabolite_col = drug_metabolite
)
