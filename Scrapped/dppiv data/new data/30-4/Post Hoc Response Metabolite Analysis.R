library(dplyr)
library(broom)
library(metadeconfoundR)
library(readr)
library(fuzzyjoin)

#Input Parameters
metadata <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
chemical_annotation <- as.data.frame(read_csv("chemical annotation.csv"))
dose_only_meta <- as.data.frame(read_csv("metadata2.csv",show_col_types = FALSE))
dose_serum_meta <- as.data.frame(read_csv("metadata3.csv",show_col_types = FALSE))
metabolite_data <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))

#Preprocessing
metadata_d <- metadata %>% filter(SampleID %in% dose_only_meta$SampleID) %>% select(SampleID,AGE,GENDER,CENTER_C,GLYCATHB)
metadata_ds <- metadata %>% filter(SampleID %in% dose_serum_meta$SampleID) %>% select(SampleID,AGE,GENDER,CENTER_C,GLYCATHB)
metabolite_d <- metabolite_data %>% filter(MGS_PID %in% metadata_d$SampleID)
metabolite_ds <- metabolite_data %>% filter(MGS_PID %in% metadata_ds$SampleID)

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
metadecounfoundR_d <- perform_metadeconfound_analysis(metadata = metadata_d,metabolite_data = metabolite_d,variables = c("GLYCATHB"),filtering_criteria = c("OK_nc","OK_sd"),output_number = 1,chemical_annotation = chemical_annotation)
metadeconfoundR_ds <- perform_metadeconfound_analysis(metadata = metadata_ds,metabolite_data = metabolite_ds,variables = c("GLYCATHB"),filtering_criteria = c("OK_nc","OK_sd"),output_number = 2,chemical_annotation = chemical_annotation)
#Linear Regression
perform_and_extract_significant_results_GLYCATHB <- function(metadata_1, metabolite_data, output_number) {
  # Combine Data for Linear Regression
  combined_data <- left_join(metabolite_data, metadata_1, by = c("MGS_PID" = "SampleID"))
  combined_data <- na.omit(combined_data)
  combined_data$MGS_PID <- NULL
  combined_data <- combined_data %>% select(AGE, GENDER, CENTER_C, GLYCATHB, everything())
  
  perform_linear_regression <- function(combined_data) {
    all_p_values <- data.frame()
    
    # Ensure column names are unique
    colnames(combined_data) <- make.unique(colnames(combined_data))
    
    for (i in 5:ncol(combined_data)) {  # Metabolites start from column 5 to the last column
      metabolite <- colnames(combined_data)[i]
      metabolite <- trimws(metabolite)  # Trim any extra spaces
      # Construct the formula correctly with quotes around the metabolite name
      formula <- as.formula(paste("GLYCATHB ~ `", metabolite, "` + AGE + GENDER + CENTER_C", sep = ""))
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
  
  # Perform linear regression for GLYCATHB
  significant_results <- perform_linear_regression(combined_data)
  
  # Export results to CSV file
  output_file <- paste0("GLYCATHB_", output_number, "_lm.csv")
  write_csv(significant_results, output_file)
  
  return(significant_results)
}

lm_d <- perform_and_extract_significant_results_GLYCATHB(metadata_1 = metadata_d,metabolite_data = metabolite_d,output_number = 1)
lm_ds <- perform_and_extract_significant_results_GLYCATHB(metadata_1 = metadata_ds,metabolite_data = metabolite_ds,output_number = 2)
