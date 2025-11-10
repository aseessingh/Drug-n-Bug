library(readr)
library(dplyr)

treated_group <- read_csv("group 3.csv" , show_col_types = FALSE)
untreated_group <- read_csv("group 2.csv" , show_col_types = FALSE)
healthy_group <- read_csv("group 8.csv",show_col_types = FALSE)
dataset <- read_csv("hub.cellcount.stats.v3.mgsamples.r.csv", show_col_types = FALSE)
drug <- c("SU_C")
#Filter for drug intake
treated_group<- treated_group %>% filter(.data[[drug]] != 0)
#Extract Sample IDs
treated_group <- treated_group %>% select(SampleID)
untreated_group <- untreated_group %>% select(SampleID)
healthy_group <- healthy_group %>% select(SampleID)
#Fix Metagenomic Sample IDs
dataset$MGSampleID <- sub("M0_","",dataset$MGSampleID)
#Adding metagenomic parameters
treated_group <- left_join(treated_group,dataset, by= c("SampleID"="MGSampleID"))
untreated_group <- left_join(untreated_group,dataset, by= c("SampleID"="MGSampleID"))
healthy_group <- left_join(healthy_group,dataset, by= c("SampleID"="MGSampleID"))
#Processing and Cleaning
treated_group <- treated_group %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)
untreated_group <- untreated_group  %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)
healthy_group <- healthy_group  %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)
treated_group <- as.data.frame(treated_group)
untreated_group <- as.data.frame(untreated_group)
healthy_group <- as.data.frame(healthy_group)
#Mean Analysis per group
output_df <- data.frame(
  Parameter = c("Shannon Diversity (Genus)", "Observed Genus Richness", "Observed Specie Richness", "Microbial Load", "Shannon Diversity (Species)", "Microbial Load Index"),
  
  # Calculate means for Treated_Group
  Treated_Group = c(mean(treated_group[,2]), mean(treated_group[,3]), mean(treated_group[,4]), mean(treated_group[,5]), mean(treated_group[,6]), mean(treated_group[,7])),
  
  # Calculate means for Untreated_Group
  Untreated_Group = c(mean(untreated_group[,2]), mean(untreated_group[,3]), mean(untreated_group[,4]), mean(untreated_group[,5]), mean(untreated_group[,6]), mean(untreated_group[,7])),
  
  # Calculate means for Healthy_Group
  Healthy_Group = c(mean(healthy_group[,2]), mean(healthy_group[,3]), mean(healthy_group[,4]), mean(healthy_group[,5]), mean(healthy_group[,6]), mean(healthy_group[,7]))
)
write_csv(output_df,"all_groups.csv")
