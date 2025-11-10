library(dplyr)
library(readr)
library(ggplot2)
library(broom) 

metadata <- read.csv("/Volumes/as13324/home/metadrugs2/MetaCardis Data Hub/hub.metadata.samples.v10.r.csv")
dose_data <- read.csv("/Volumes/as13324/home/metadrugs2/MetaCardis Data Hub/MC_Drug_doses.csv")
metabolite_data <- read_csv("/Volumes/as13324/home/metadrugs2/MetaCardis Data Hub/MC_RawMetabolonBatch_norm_ImputedDatawoDups_KC_141223.csv")

drug <- "METFORMIN_C"
drug_metabolite <- "metformin"
dose_drug <- "DOSAGE_METFORMIN_C"

# Preprocess metadata
metadata <- metadata %>% 
  filter(!!sym(drug) == 1) %>%
  select(SampleID, AGE, GENDER, CENTER_C, BMI_C)

metadata$GENDER <- as.factor(metadata$GENDER)
metadata$CENTER_C <- as.factor(metadata$CENTER_C)

# Select metabolite data
metabolite_data <- metabolite_data %>% select(all_of(drug_metabolite), MGS_PID)

# Select Drug Dosing Data
dose_data <- dose_data %>% select(all_of(dose_drug), SampleID)

# Merge metadata and metabolite data
merged_metadata_metabolite <- inner_join(metadata, metabolite_data, by = c("SampleID" = "MGS_PID"))
merged_dose_meta_metabolite <- inner_join(merged_metadata_metabolite, dose_data, by = "SampleID")

# Step 1: Remove rows with non-finite values in the response column
merged_dose_meta_metabolite <- merged_dose_meta_metabolite %>%
  filter(!is.na(!!sym(drug_metabolite)) & !is.nan(!!sym(drug_metabolite)) & !is.infinite(!!sym(drug_metabolite)))

# Step 2: Replace string "NA" with actual NA values, only for character columns
merged_dose_meta_metabolite <- merged_dose_meta_metabolite %>%
  mutate(across(where(is.character), ~ na_if(., "NA")))

# Step 3: Filter out rows where more than 3 columns have NA values
merged_dose_meta_metabolite <- merged_dose_meta_metabolite %>%
  filter(rowSums(is.na(.)) <= 1)

# Fit the models
formula0 <- reformulate(dose_drug, response = drug_metabolite)
model0 <- lm(formula = formula0, data = merged_dose_meta_metabolite)

formula1 <- reformulate(c(dose_drug, "AGE", "GENDER", "CENTER_C"), response = drug_metabolite)
model1 <- lm(formula = formula1, data = merged_dose_meta_metabolite)

formula2 <- reformulate(c(dose_drug, "AGE", "GENDER", "CENTER_C", "BMI_C"), response = drug_metabolite)
model2 <- lm(formula = formula2, data = merged_dose_meta_metabolite)

# Create a dataframe to store model summaries
model_summaries <- list(summary(model0), summary(model1), summary(model2))

tidy_summaries <- lapply(model_summaries, tidy)
all_model_summaries <- bind_rows(tidy_summaries, .id = "model_id")

write.csv(all_model_summaries, "model_summaries.csv", row.names = FALSE)

# Function to create plots
plot_model <- function(model, data, model_id) {
  predictions <- predict(model, newdata = data)
  ggplot(data, aes(x = !!sym(drug_metabolite), y = predictions)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    labs(title = paste("Model", model_id, "Prediction vs Actual"),
         x = paste("Actual", drug_metabolite),
         y = paste("Predicted", drug_metabolite)) +
    theme_minimal()
}

plot0 <- plot_model(model0, merged_dose_meta_metabolite, 0)
plot1 <- plot_model(model1, merged_dose_meta_metabolite, 1)
plot2 <- plot_model(model2, merged_dose_meta_metabolite, 2)

ggsave("model0_plot.png", plot0)
ggsave("model1_plot.png", plot1)
ggsave("model2_plot.png", plot2)

# Predictions
predictions_model0 <- predict(model0, newdata = merged_dose_meta_metabolite)
predictions_model1 <- predict(model1, newdata = merged_dose_meta_metabolite)
predictions_model2 <- predict(model2, newdata = merged_dose_meta_metabolite)

# Consolidated plot
consolidated_plot <- ggplot() +
  geom_point(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model0), color = "blue") + 
  geom_smooth(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model0), method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_point(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model1), color = "red") + 
  geom_smooth(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_point(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model2), color = "green") + 
  geom_smooth(data = merged_dose_meta_metabolite, aes(x = !!sym(drug_metabolite), y = predictions_model2), method = "lm", se = FALSE, color = "green", linetype = "dashed") +
  labs(title = "Comparison of Predictions from All Models",
       x = paste("Actual", drug_metabolite),
       y = paste("Predicted", drug_metabolite)) +
  theme_minimal()

ggsave("consolidated_plot.png", consolidated_plot)
