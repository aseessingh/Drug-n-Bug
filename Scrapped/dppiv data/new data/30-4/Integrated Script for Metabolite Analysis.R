library(dplyr)
library(readr)
library(metadeconfoundR)
library(fuzzyjoin)
library(broom)
library(tibble)

#Input for Code
metabolite_data <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))
metabolite_raw_for_drug <- as.data.frame(read_csv("MC_raw_metabolites.csv",show_col_types = FALSE))
dose_data <- as.data.frame(read_csv("MC_Drug_doses.csv",show_col_types = FALSE))
metadata <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
MGS <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv" , show_col_types = FALSE))
drug_dose <- c("DOSAGE_DPPIV_C")
drug_metabolite <- c("sitagliptin")
drug_metadata <- c("DPPIV_C")
chemical_annotation <- as.data.frame(read_csv("chemical annotation.csv"))

#Preprocessing to get first 3 datasets
MGS$MGSampleID <- gsub("M0_","",MGS$MGSampleID)
metadata <- metadata %>% select(SampleID, AGE, GENDER, CENTER_C,drug_metadata)
metadata1_ids <- intersect(metadata$SampleID,metabolite_data$MGS_PID)
metadata1_ids <- intersect(metadata1_ids,MGS$MGSampleID)
metadata_1 <- metadata %>% filter(SampleID %in% metadata1_ids,get(drug_metadata) != 0) %>% select(SampleID,AGE,GENDER,CENTER_C)
metabolite_data_1 <- metabolite_data %>% filter(MGS_PID %in% metadata1_ids)
drug_dose_df <- dose_data %>% filter(get(drug_dose) !=0) %>% select(MCid, drug_dose)
serum_drug_df <- metabolite_raw_for_drug %>% filter(get(drug_metabolite) !=0) %>% select(MGS_PID,drug_metabolite)
metadata_1 <- left_join(metadata_1,drug_dose_df,by=c("SampleID"="MCid"))
metadata_1 <- left_join(metadata_1,serum_drug_df,by=c("SampleID"="MGS_PID"))
metadata2_ids <- intersect(metadata1_ids,drug_dose_df$MCid)
metadata_2 <- metadata_1 %>% filter(SampleID %in% metadata2_ids)
metabolite_data_2 <- metabolite_data_1 %>% filter(MGS_PID %in% metadata2_ids)
metadata3_ids <- intersect(metadata2_ids,serum_drug_df$MGS_PID)
metadata_3 <- metadata_2 %>% filter(SampleID %in% metadata3_ids)
metabolite_data_3 <- metabolite_data_2 %>% filter(MGS_PID %in% metadata3_ids)
metadata_4 <- metadata %>% filter(get(drug_metadata)==0,SampleID %in% MGS$MGSampleID) %>% select(SampleID, AGE, GENDER, CENTER_C)
serum_drug_anti_df <- metabolite_raw_for_drug %>% filter(get(drug_metabolite)==0) %>% select(MGS_PID,drug_metabolite)
dose_drug_anti_df <- dose_data %>% filter(get(drug_dose)==0) %>% select(MCid,drug_dose)
metadata_4 <- left_join(metadata_4,serum_drug_anti_df,by=c("SampleID"="MGS_PID"))
metadata_4 <- left_join(metadata_4,dose_drug_anti_df,by=c("SampleID"="MCid"))
metadata_4_lm <- metadata_4
metadata_3_D <- metadata_3
metadata_3_D$Status <- 1
metadata_4$Status <- 0
metadata_3_D <- metadata_3_D %>% select(Status,everything())
metadata_4 <- metadata_4 %>% select(SampleID,Status,everything())
metadata_4 <- rbind(metadata_4,metadata_3_D)
metadata_4 <- metadata_4 %>% mutate(across(-1, ~ ifelse(is.na(.x), 0, .x)))
metabolite_data_4 <- metabolite_data %>% filter(MGS_PID %in% metadata_4$SampleID)
write_csv(metadata_1,"metadata1.csv")
write_csv(metadata_2,"metadata2.csv")
write_csv(metadata_3,"metadata3.csv")
write_csv(metadata_4,"metadata4.csv")

#MetadeconfoundR
perform_metadeconfound_analysis <- function(metadata, metabolite_data, variables, filtering_criteria, chemical_annotation, output_number, additional_continuous = NULL, additional_categorical = NULL) {
  # Default continuous and categorical variables
  type_continuous <- c("AGE")
  type_categorical <- c("GENDER", "CENTER_C")
  
  # Add additional continuous and categorical variables if provided
  if (!is.null(additional_continuous)) {
    type_continuous <- c(type_continuous, additional_continuous)
  }
  if (!is.null(additional_categorical)) {
    type_categorical <- c(type_categorical, additional_categorical)
  }
  
  # Preprocessing
  rownames(metadata) <- metadata[, 1]
  rownames(metabolite_data) <- metabolite_data[, 1]
  metadata[, 1] <- NULL
  metabolite_data[, 1] <- NULL
  metadata <- metadata[match(rownames(metabolite_data), rownames(metadata)), ]
  metabolite_data <- metabolite_data[match(rownames(metadata), rownames(metabolite_data)), ]
  
  # Check if "Status" is in additional categorical variables
  contains_status <- any(grepl("Status", type_categorical))
  
  # Naive MetaDeconfound analysis
  naive_result <- MetaDeconfound(
    featureMat = metabolite_data,
    metaMat = metadata,
    QCutoff = 0.05,
    DCutoff = 0.2,
    typeContinuous = type_continuous,
    typeCategorical = type_categorical,
    robustCutoffRho = 1,
    startStop = "naive",
    adjustMethod = "fdr",
    robustCutoff = ifelse(contains_status, NULL, 0)
  )
  
  # Posthoc MetaDeconfound analysis
  posthoc_result <- MetaDeconfound(
    featureMat = metabolite_data,
    metaMat = metadata,
    QValues = naive_result$Qs,
    DValues = naive_result$Ds,
    QCutoff = 0.05,
    DCutoff = 0.2,
    typeContinuous = type_continuous,
    typeCategorical = type_categorical,
    robustCutoffRho = 1,
    robustCutoff = ifelse(contains_status, NULL, 0)
  )
  
  # Extract and filter results
  result <- as.data.frame(cbind(naive_result$Ps, naive_result$Qs, naive_result$Ds, posthoc_result$status))
  result$FuzzyName <- rownames(result)
  colnames(result) <- make.unique(colnames(result))
  
  # Use the chemical_annotation data frame directly
  max_distance <- 20
  method <- "osa"
  
  # Annotate results
  annotate_until_complete <- function(adj, chemical_annotation, max_distance = 20, method = "osa") {
    annotated_all <- tibble()
    remaining <- adj
    dist <- 0
    
    while (nrow(remaining) > 0 && dist <= max_distance) {
      message("Annotating with max_dist = ", dist)
      
      annotated_current <- stringdist_inner_join(
        remaining,
        chemical_annotation,
        by = c("FuzzyName" = "CHEMICAL_NAME"),
        method = method,
        max_dist = dist
      )
      
      if (nrow(annotated_current) == 0) {
        dist <- dist + 1
        next
      }
      
      annotated_all <- bind_rows(annotated_all, annotated_current) %>%
        distinct(`FuzzyName`, .keep_all = TRUE)
      
      remaining <- adj %>% filter(!(`FuzzyName` %in% annotated_all$`FuzzyName`))
      dist <- dist + 1
    }
    
    return(annotated_all)
  }
  
  result <- annotate_until_complete(adj = result, chemical_annotation, max_distance, method)
  
  # Iterate over variables and filter results
  for (variable in variables) {
    column_name <- paste0(variable, ".3")
    filtered_result <- result %>% filter(get(column_name) %in% filtering_criteria)
    
    # Automatically name the output file based on the variable and output number
    output_file <- paste0(variable, output_number, ".csv")
    
    # Ensure the output file path is a valid string
    if (is.character(output_file) && nzchar(output_file)) {
      write_csv(filtered_result, output_file)
    } else {
      stop("Invalid output file path: ", output_file)
    }
  }
  
  return(result)
}
metadeconfound_1 <- perform_metadeconfound_analysis(metadata = metadata_1,metabolite_data = metabolite_data_1,variable = c(drug_dose,drug_metabolite),chemical_annotation = chemical_annotation,filtering_criteria = c("OK_nc","OK_sd"),output_number = 1)
metadeconfound_2 <- perform_metadeconfound_analysis(metadata = metadata_2,metabolite_data = metabolite_data_2,variable = c(drug_dose,drug_metabolite),chemical_annotation = chemical_annotation,filtering_criteria = c("OK_nc","OK_sd"),output_number = 2)
metadeconfound_3 <- perform_metadeconfound_analysis(metadata = metadata_3,metabolite_data = metabolite_data_3,variable = c(drug_dose,drug_metabolite),chemical_annotation = chemical_annotation,filtering_criteria = c("OK_nc","OK_sd"),output_number = 3)
metadeconfound_4 <- perform_metadeconfound_analysis(metadata = metadata_4,metabolite_data = metabolite_data_4,variable = c(drug_dose,drug_metabolite),chemical_annotation = chemical_annotation,filtering_criteria = c("OK_nc","OK_sd"),output_number = 4)

#Linear Regression
perform_and_extract_significant_results <- function(metadata_1, metabolite_data, drug_dose, drug_metabolite, output_number) {
  # Split metadata into dose and serum
  metadata_dose <- metadata_1 %>% select(SampleID, AGE, GENDER, CENTER_C, drug_dose)
  metadata_serum <- metadata_1 %>% select(SampleID, AGE, GENDER, CENTER_C, drug_metabolite)
  
  # Combine Data for Linear Regression
  combined_data_dose <- left_join(metabolite_data, metadata_dose, by = c("MGS_PID" = "SampleID"))
  combined_data_dose <- na.omit(combined_data_dose)
  combined_data_dose$MGS_PID <- NULL
  combined_data_dose <- combined_data_dose %>% select(AGE, GENDER, CENTER_C, drug_dose, everything())
  
  combined_data_serum <- left_join(metabolite_data, metadata_serum, by = c("MGS_PID" = "SampleID"))
  combined_data_serum <- na.omit(combined_data_serum)
  combined_data_serum$MGS_PID <- NULL
  combined_data_serum <- combined_data_serum %>% select(AGE, GENDER, CENTER_C, drug_metabolite, everything())
  
  perform_linear_regression <- function(combined_data, dependent_var) {
    all_p_values <- data.frame()
    
    # Ensure column names are unique
    colnames(combined_data) <- make.unique(colnames(combined_data))
    
    for (i in 5:1362) {  # Metabolites start from column 5 to 1362
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
    all_p_values$adjusted_p_value <- p.adjust(all_p_values$p_value, method = "fdr")
    
    # Extract significant results
    significant_results <- all_p_values %>% filter(adjusted_p_value < 0.05)
    
    return(significant_results)
  }
  
  # Perform linear regression for dose and serum data
  significant_results_dose <- perform_linear_regression(combined_data_dose, drug_dose)
  significant_results_serum <- perform_linear_regression(combined_data_serum, drug_metabolite)
  
  # Export results to CSV files
  output_file_dose <- paste0(drug_dose, output_number, "_lm.csv")
  output_file_serum <- paste0(drug_metabolite, output_number, "_lm.csv")
  
  write_csv(significant_results_dose, output_file_dose)
  write_csv(significant_results_serum, output_file_serum)
  
  return(list(dose_data_results = significant_results_dose, drug_results_serum = significant_results_serum))
}

lm_1 <- perform_and_extract_significant_results(metadata_1 = metadata_1,metabolite_data = metabolite_data_1,drug_dose = drug_dose,drug_metabolite = drug_metabolite,output_number = 1)
lm_2 <- perform_and_extract_significant_results(metadata_1 = metadata_2,metabolite_data = metabolite_data_2,drug_dose = drug_dose,drug_metabolite = drug_metabolite,output_number = 2)
lm_3 <- perform_and_extract_significant_results(metadata_1 = metadata_3,metabolite_data = metabolite_data_3,drug_dose = drug_dose,drug_metabolite = drug_metabolite,output_number = 3)
lm_4 <- perform_and_extract_significant_results(metadata_1 = metadata_4,metabolite_data = metabolite_data_4,drug_dose = drug_dose,drug_metabolite = drug_metabolite,output_number = 4)

