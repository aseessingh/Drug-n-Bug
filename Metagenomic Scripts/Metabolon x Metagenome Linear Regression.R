library(dplyr)
library(readr)

#Input
metadata <- as.data.frame(read_csv("metadata3.csv",show_col_types = FALSE))
group_3 <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
metabolite <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))
metagenome <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv",show_col_types = FALSE))
results_metagenome <- as.data.frame(read_csv("metformin vs non-metformin MGS.csv",show_col_types = FALSE))
results_metabolon <- as.data.frame(read_csv("lm_AG.csv",show_col_types = FALSE))
dose_data <- c("DOSAGE_METFORMIN_C")
drug_metabolite <- c("metformin")

#Preprocessing
group_3 <- group_3 %>% select(SampleID,GLYCATHB)
metadata <- left_join(metadata,group_3,by=c("SampleID"))
metabolite_for_AG <- metabolite
metabolite_for_AG <- metabolite_for_AG %>% select(MGS_PID,`1,5-anhydroglucitol (1,5-AG)`) %>% rename (AG=`1,5-anhydroglucitol (1,5-AG)`)
metadata <- left_join(metadata,metabolite_for_AG,by=c("SampleID"="MGS_PID"))
metagenome <- metagenome %>% select(MGSampleID,results_metagenome$Feature)
metagenome$MGSampleID <- gsub("M0_","",metagenome$MGSampleID)
metadata <- left_join(metadata,metagenome,by=c("SampleID"="MGSampleID"))
metabolites_AG <- metabolite %>% select(MGS_PID,results_metabolon$metabolite)
metadata <- left_join(metadata,metabolites_AG,by=c("SampleID"="MGS_PID"))

#Lm
lm_analysis_both_ways <- function(df, metagenome_cols, metabolon_cols) {
  library(dplyr)
  library(broom)
  library(tibble)
  
  # Function to adjust p-values using FDR
  adjust_pvalues <- function(pvalues) {
    p.adjust(pvalues, method = "fdr")
  }
  
  perform_linear_regression <- function(combined_data, dependent_var, independent_vars) {
    all_p_values <- data.frame()
    
    # Ensure column names are unique
    colnames(combined_data) <- make.unique(colnames(combined_data))
    
    for (i in independent_vars) {
      independent_var <- colnames(combined_data)[i]
      independent_var <- trimws(independent_var)  # Trim any extra spaces
      # Construct the formula correctly with quotes around the variable names
      formula <- as.formula(paste("`", dependent_var, "` ~ `", independent_var, "` + AGE + GENDER + CENTER_C", sep = ""))
      model <- lm(formula, data = combined_data)
      summary_model <- summary(model)
      
      p_value <- summary_model$coefficients[2, "Pr(>|t|)"]  # The independent variable is always the second coefficient
      beta_coefficient <- summary_model$coefficients[2, "Estimate"]
      all_p_values <- rbind(all_p_values, data.frame(dependent_var = dependent_var, independent_var = independent_var, p_value = p_value, beta_coefficient = beta_coefficient))
    }
    
    # Apply FDR correction
    all_p_values$adjusted_p_value <- adjust_pvalues(all_p_values$p_value)
    
    # Extract significant results
    significant_results <- all_p_values %>% filter(adjusted_p_value < 0.05)
    
    return(significant_results)
  }
  
  # Perform linear regression for metagenome ~ metabolon + AGE + GENDER + CENTER_C
  results_metagenome <- data.frame()
  for (metagenome_col in metagenome_cols) {
    results_metagenome <- rbind(results_metagenome, perform_linear_regression(df, colnames(df)[metagenome_col], metabolon_cols))
  }
  
  # Perform linear regression for metabolon ~ metagenome + AGE + GENDER + CENTER_C
  results_metabolon <- data.frame()
  for (metabolon_col in metabolon_cols) {
    results_metabolon <- rbind(results_metabolon, perform_linear_regression(df, colnames(df)[metabolon_col], metagenome_cols))
  }
  
  # Export results to CSV files
  write.csv(results_metagenome, "Metabolon_x_Metagenome.csv", row.names = FALSE)
  write.csv(results_metabolon, "Metagenome_x_Metabolon.csv", row.names = FALSE)
  
  return(list(
    metagenome_results = results_metagenome,
    metabolon_results = results_metabolon
  ))
}
lm <- lm_analysis_both_ways(df=metadata, metagenome_cols = 9:18,metabolon_cols = 19:349)
