library(dplyr)
library(readr)
library(ggplot2)
library(broom) 

metadata <- hub_metadata_samples_v10_r 
metabolite_data <- MC_RawMetabolonBatch_norm_ImputedDatawoDups_KC_141223
drug <- c("METFORMIN_C")
drug_metabolite <- c("metformin")
response <- c("GLYCATHB")

# Preprocess metadata
metadata <- hub_metadata_samples_v10_r %>% 
  filter(!!sym(drug) == 1) %>%
  select(SampleID,AGE, GENDER, CENTER_C, response, BMI_C)

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
  filter(rowSums(is.na(.)) <= 1)

# Fit the models
formula0 <- reformulate(drug_metabolite, response = response)  # Replaced GLYCATHB with 'response'
model0 <- lm(formula=formula0, data=merged_metadata_metabolite)

formula1 <- reformulate(termlabels = c(drug_metabolite, "AGE", "GENDER", "CENTER_C"), response = response)  # Replaced GLYCATHB with 'response'
model1 <- lm(formula=formula1, data=merged_metadata_metabolite)

formula2 <- reformulate(termlabels = c(drug_metabolite, "AGE", "GENDER", "CENTER_C", "BMI_C"), response = response)  # Replaced GLYCATHB with 'response'
model2 <- lm(formula=formula2, data=merged_metadata_metabolite)

# Create a dataframe to store model summaries
model_summaries <- list(
  summary(model0),
  summary(model1),
  summary(model2)
)

# Extract tidy summaries from each model
tidy_summaries <- lapply(model_summaries, tidy)

# Combine all summaries into one data frame
all_model_summaries <- bind_rows(tidy_summaries, .id = "model_id")

# Export summaries to CSV
write.csv(all_model_summaries, "model_summaries.csv", row.names = FALSE)

# Function to create plots
plot_model <- function(model, data, model_id) {
  predictions <- predict(model, newdata = data)
  ggplot(data, aes(x = !!sym(response), y = predictions)) +  # Replaced GLYCATHB with 'response'
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    labs(title = paste("Model", model_id, "Prediction vs Actual"),
         x = paste("Actual", response),
         y = paste("Predicted", response)) +
    theme_minimal()
}

# Create individual plots for each model
plot0 <- plot_model(model0, merged_metadata_metabolite, 0)
plot1 <- plot_model(model1, merged_metadata_metabolite, 1)
plot2 <- plot_model(model2, merged_metadata_metabolite, 2)

# Save the individual plots
ggsave("model0_plot.png", plot0)
ggsave("model1_plot.png", plot1)
ggsave("model2_plot.png", plot2)

# Generate predictions for all models using the same dataset
predictions_model0 <- predict(model0, newdata = merged_metadata_metabolite)
predictions_model1 <- predict(model1, newdata = merged_metadata_metabolite)
predictions_model2 <- predict(model2, newdata = merged_metadata_metabolite)

# Now, plot the results
consolidated_plot <- ggplot() +
  geom_point(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model0), color = "blue") + 
  geom_smooth(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model0), method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_point(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model1), color = "red") + 
  geom_smooth(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_point(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model2), color = "green") + 
  geom_smooth(data = merged_metadata_metabolite, aes(x = !!sym(response), y = predictions_model2), method = "lm", se = FALSE, color = "green", linetype = "dashed") +
  labs(title = "Comparison of Predictions from All Models",
       x = paste("Actual", response),
       y = paste("Predicted", response)) +
  theme_minimal()

# Save the consolidated plot
ggsave("consolidated_plot.png", consolidated_plot)
