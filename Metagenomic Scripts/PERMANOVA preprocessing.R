library(readr)
library(dplyr)
library(DescTools)

treated_group <- read_csv("group 3.csv" , show_col_types = FALSE)
untreated_group <- read_csv("group 2.csv" , show_col_types = FALSE)
healthy_group <- read_csv("group 8.csv",show_col_types = FALSE)
dataset <- read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv", show_col_types = FALSE)
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
write_csv(healthy_group, "HC.csv")
write_csv(untreated_group, "UMMC.csv")
write_csv(treated_group_drug1,"TG_MET.csv")
write_csv(treated_group_drug2, "TG_DPP.csv")
write_csv(treated_group_drug3, "TG_SU.csv")
