library(dplyr)
library(readr)
library(ggplot2)
library(broom) 

metadata <- read_csv("hub.metadata.samples.v10.r.csv",show_col_types=FALSE)
metabolite_data <- read_csv("MC_RawMetabolonBatch_norm_ImputedDatawoDups_KC_141223.csv",show_col_types = FALSE)
drug <- c("METFORMIN_C")
drug_metabolite <- c("metformin")
response <- c("GLYCATHB")

# Preprocess metadata
metadata <- metadata %>% 
  filter(!!sym(drug) == 1) %>%
  select(SampleID,AGE, GENDER, CENTER_C, all_of(response), BMI_C)

metadata$GENDER <- as.factor(metadata$GENDER)
metadata$CENTER_C <- as.factor(metadata$CENTER_C)

# Select metabolite data
metabolite_data <- metabolite_data %>% select(all_of(drug_metabolite),MGS_PID)

# Merge metadata and metabolite data
merged_metadata_metabolite <- inner_join(metadata, metabolite_data, by = c("SampleID" = "MGS_PID"))

# Step 1: Remove rows with non-finite values in the 'response' column
merged_metadata_metabolite <- merged_metadata_metabolite %>%
  filter(!is.na(!!sym(response)) & !is.nan(!!sym(response)) & !is.infinite(!!sym(response)))

# Step 2: Replace string "NA" with actual NA values, only for character columns
merged_metadata_metabolite <- merged_metadata_metabolite %>%
  mutate(across(where(is.character), ~ na_if(., "NA")))  # Replace "NA" string with NA only in character columns

# Step 3: Filter out rows where more than 3 columns have NA values
merged_metadata_metabolite <- merged_metadata_metabolite %>%
  filter(rowSums(is.na(.)) <= 0)

formula2 <- reformulate(termlabels = c(drug_metabolite, "AGE", "GENDER", "CENTER_C", "BMI_C"), response = response)  # Replaced GLYCATHB with 'response'
model2 <- lm(formula=formula2, data=merged_metadata_metabolite)
model2_with_residuals <- augment(model2,data=merged_metadata_metabolite)
resid_sd <- sd(model2_with_residuals$.resid)
outlier_indices <- which(abs(model2_with_residuals$.resid) > 2 * resid_sd)
outlier_sample_ids <- model2_with_residuals[outlier_indices, "SampleID"]
write.csv(outlier_sample_ids, "outlier_sample_ids.csv", row.names = FALSE)

