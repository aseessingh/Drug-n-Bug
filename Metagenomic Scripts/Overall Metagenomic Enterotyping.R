library(readr)
library(dplyr)
library(DescTools)

treated_group <- read_csv("group 3.csv" , show_col_types = FALSE)
untreated_group <- read_csv("group 2.csv" , show_col_types = FALSE)
healthy_group <- read_csv("group 8.csv",show_col_types = FALSE)
dataset <- read_csv("hub.enterotypeHellinger.mgsamples.v4.r.csv", show_col_types = FALSE)
drug <- c("METFORMIN_C", "DPPIV_C", "SU_C")
#Untreated and Healthy Group Data
untreated_group <- untreated_group %>% select(SampleID)
healthy_group <- healthy_group %>% select(SampleID)
#Fix Metagenomic Sample IDs
dataset$MGSampleID <- sub("M0_","",dataset$MGSampleID)
#Adding metagenomic parameters
treated_group_drug1 <- as.data.frame(
  treated_group %>%
    filter(.data[[drug[1]]] != 0) %>%
    select(SampleID) %>%
    left_join(dataset, by = c("SampleID" = "MGSampleID")) %>%
    filter(rowSums(is.na(.) | . == "NA") == 0) %>%
    distinct(SampleID, .keep_all = TRUE)
)

treated_group_drug2 <- as.data.frame(
  treated_group %>%
    filter(.data[[drug[2]]] != 0) %>%
    select(SampleID) %>%
    left_join(dataset, by = c("SampleID" = "MGSampleID")) %>%
    filter(rowSums(is.na(.) | . == "NA") == 0) %>%
    distinct(SampleID, .keep_all = TRUE)
)

treated_group_drug3 <- as.data.frame(
  treated_group %>%
    filter(.data[[drug[3]]] != 0) %>%
    select(SampleID) %>%
    left_join(dataset, by = c("SampleID" = "MGSampleID")) %>%
    filter(rowSums(is.na(.) | . == "NA") == 0) %>%
    distinct(SampleID, .keep_all = TRUE)
)

untreated_group <- left_join(untreated_group,dataset, by= c("SampleID"="MGSampleID"))
healthy_group <- left_join(healthy_group,dataset, by= c("SampleID"="MGSampleID"))
#Processing and Cleaning
untreated_group <- untreated_group  %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)
healthy_group <- healthy_group  %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)
untreated_group <- as.data.frame(untreated_group)
healthy_group <- as.data.frame(healthy_group)
#Mean Analysis per group
output_df <- data.frame(
  Parameter = c("Most Common Enterotype"),
  Healthy_Group = c(Mode(healthy_group[,2])),
  Untreated_Group = c(Mode(untreated_group[,2])),
  Drug1 = c(Mode(treated_group_drug1[,2])),
  Drug2 = c(Mode(treated_group_drug2[,2])),
  Drug3 = c(Mode(treated_group_drug3[,2]))
)
colnames(output_df)[4:6] <- drug
write_csv(output_df,"all_groups.csv")
