library(readr)
library(dplyr)
library(purrr)
merge_csvs_metad <- function(file1, file2, file3, drug_dose, drug_metabolite) {
  # Read the CSV files
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  df3 <- read.csv(file3)
  
  # Define the columns with .2 suffix
  df1$drug_dose.2 <- df1[[drug_dose]]
  df2$drug_metabolite.2 <- df2[[drug_metabolite]]
  
  # Select and rename the required columns
  df1 <- df1 %>%
    select(drug_dose.2, CHEMICAL_NAME, ontology) %>%
    rename(Corr_Dose = drug_dose.2)
  
  df2 <- df2 %>%
    select(drug_metabolite.2, CHEMICAL_NAME, ontology) %>%
    rename(Corr_Serum = drug_metabolite.2)
  
  df3 <- df3 %>%
    select(GLYCATHB, CHEMICAL_NAME, ontology) %>%
    rename(Corr_HbA1c = GLYCATHB)
  
  # Merge the dataframes
  merged_df <- reduce(list(df1, df2, df3), function(x, y) merge(x, y, by = c("CHEMICAL_NAME", "ontology"), all = TRUE))
  
  return(merged_df)
}

metad_1 <- merge_csvs_metad("DOSAGE_METFORMIN_C1.csv","metformin1.csv","GLYCATHB1.csv","DOSAGE_METFORMIN_C","metformin")
metad_2 <- merge_csvs_metad("DOSAGE_METFORMIN_C2.csv","metformin2.csv","GLYCATHB1.csv","DOSAGE_METFORMIN_C","metformin")
metad_3 <- merge_csvs_metad("DOSAGE_METFORMIN_C3.csv","metformin3.csv","GLYCATHB2.csv","DOSAGE_METFORMIN_C","metformin")
metad_4 <- merge_csvs_metad("DOSAGE_METFORMIN_C4.csv","metformin4.csv","GLYCATHB2.csv","DOSAGE_METFORMIN_C","metformin")
write.csv(metad_1,"metad_case1.csv")
write.csv(metad_2,"metad_case2.csv")
write.csv(metad_3,"metad_case3.csv")
write.csv(metad_4,"metad_case4.csv")

merge_csvs_lm <- function(file1, file2, file3) {
  # Helper function to extract the new column name from the file name
  get_new_col_name <- function(file) {
    base_name <- tools::file_path_sans_ext(basename(file))
    if (grepl("DOSAGE_DPPIV", base_name)) {
      return("Beta_Dose")
    } else if (grepl("sitagliptin", base_name)) {
      return("Beta_Serum")
    } else if (grepl("glycathb", base_name, ignore.case = TRUE)) {
      return("Beta_HbA1c")
    } else {
      stop("Unknown file type")
    }
  }
  
  # Read the CSV files
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  df3 <- read.csv(file3)
  
  # Extract new column names
  col_name1 <- get_new_col_name(file1)
  col_name2 <- get_new_col_name(file2)
  col_name3 <- get_new_col_name(file3)
  
  # Select and rename the required columns
  df1 <- df1 %>%
    select(metabolite, beta_coefficient) %>%
    rename(!!col_name1 := beta_coefficient)
  
  df2 <- df2 %>%
    select(metabolite, beta_coefficient) %>%
    rename(!!col_name2 := beta_coefficient)
  
  df3 <- df3 %>%
    select(metabolite, beta_coefficient) %>%
    rename(!!col_name3 := beta_coefficient)
  
  # Merge the dataframes
  merged_df <- reduce(list(df1, df2, df3), function(x, y) merge(x, y, by = "metabolite", all = TRUE))
  
  return(merged_df)
}

lm_1 <- merge_csvs_lm("DOSAGE_DPPIV_C1_lm.csv","sitagliptin1_lm.csv","GLYCATHB_1_lm.csv")
lm_2 <- merge_csvs_lm("DOSAGE_DPPIV_C2_lm.csv","sitagliptin2_lm.csv","GLYCATHB_1_lm.csv")
lm_3 <- merge_csvs_lm("DOSAGE_DPPIV_C3_lm.csv","sitagliptin3_lm.csv","GLYCATHB_2_lm.csv")
lm_4 <- merge_csvs_lm("DOSAGE_DPPIV_C4_lm.csv","sitagliptin4_lm.csv","GLYCATHB_2_lm.csv")
write.csv(lm_1,"lm_case1.csv")
write.csv(lm_2,"lm_Case2.csv")
write.csv(lm_3,"lm_case3.csv")
write.csv(lm_4,"lm_case4.csv")

