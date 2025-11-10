library(dplyr)
library(metadeconfoundR)
library(readr)

#Input
metabolite_data_without_drugs <- read_csv("MC_metabolites_wo_Drugs.csv", show_col_types = FALSE)
treated_group <- read_csv("group 3.csv", show_col_types = FALSE)
healthy_group <- read_csv("group 8.csv", show_col_types = FALSE)
response <- c("GLYCATHB")
drug_info <- c("METFORMIN_C")

#Make Code Fit
group3 <- treated_group
group8 <- healthy_group
drug_info_metadata <- drug_info

#Filtering for Drug and Cleaning
drug_group3 <- group3 %>% filter(drug_info_metadata != 0) %>% select(all_of(c("SampleID","AGE", "GENDER", "CENTER_C",response)))
drug_metadata <- drug_group3
drug_metadata <- drug_metadata %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)

#Extract matching data from metabolite dataset
drug_metabolite <- metabolite_data_without_drugs %>% filter(MGS_PID %in% drug_metadata[[1]])
drug_metabolite <- drug_metabolite %>% distinct(MGS_PID, .keep_all = TRUE)

#Add Group 8 to datasets for control
group8_meta <- group8 %>% select(all_of(c("SampleID","AGE", "GENDER", "CENTER_C",response)))
group8_metabolites <- metabolite_data_without_drugs %>% filter(MGS_PID %in% group8_meta$SampleID)
drug_metadata$Status <- 1
group8_meta$Status <- 0
drug_metadata <- rbind(group8_meta, drug_metadata)
drug_metabolite <- rbind(drug_metabolite, group8_metabolites)
drug_metadata <- drug_metadata %>% filter(rowSums(is.na(.) | . == "NA") == 0)
#Giving Row Names and harmonizing number in both datasets
drug_metadata <- as.data.frame(drug_metadata)
drug_metabolite <- as.data.frame(drug_metabolite)
drug_metadata <- drug_metadata %>% distinct(drug_metadata$SampleID, .keep_all=TRUE)
drug_metabolite <- drug_metabolite %>% distinct(drug_metabolite$MGS_PID, .keep_all = TRUE)
drug_metadata$`drug_metadata$SampleID` <- NULL
drug_metabolite$`drug_metabolite$MGS_PID` <- NULL
drug_metadata <- drug_metadata %>% filter(SampleID %in% drug_metabolite$MGS_PID)
drug_metabolite <- drug_metabolite %>% filter(MGS_PID %in% drug_metadata$SampleID)
rownames(drug_metadata) <- drug_metadata[, 1]
rownames(drug_metabolite) <- drug_metabolite[, 1]
drug_metadata[, 1] <- NULL
drug_metabolite[, 1] <- NULL
drug_metabolite <- drug_metabolite[match(rownames(drug_metadata),rownames(drug_metabolite)),,drop=FALSE]
drug_metadata <- drug_metadata %>% select(Status, AGE, GENDER, CENTER_C, response)
write.csv(drug_metadata, "drug_metadata.csv",row.names = TRUE)
write.csv(drug_metabolite, "drug_metabolite.csv",row.names = TRUE)

#Metadeconfound Unadjusted for Demographics
drug_metadata_undaj <- drug_metadata %>% select(Status,response)
metadeconfound_naive <- MetaDeconfound(featureMat = drug_metabolite, metaMat = drug_metadata_undaj, adjustMethod = "fdr", startStop = "naiveStop", typeContinuous = response)
metadeconfound_posthoc <- MetaDeconfound(featureMat = drug_metabolite, metaMat = drug_metadata_undaj, QValues = metadeconfound_naive$Qs, DValues = metadeconfound_naive$Ds, typeContinuous = response)
unadj_metadeconfound_df <- as.data.frame(cbind(metadeconfound_naive$Ps,metadeconfound_naive$Qs,metadeconfound_naive$Ds,metadeconfound_posthoc$status))
write.csv(unadj_metadeconfound_df, "unadj_metadeconfound.csv", row.names=TRUE)
save(unadj_metadeconfound_df, file= "unadj_metadeconfound.Rdata")

#Metadeconfound Adjusted for Demographics
metadeconfound_naive_adj <- MetaDeconfound(
  featureMat = drug_metabolite,
  metaMat = drug_metadata,
  adjustMethod = "fdr",
  startStop = "naiveStop",
  typeContinuous = c(response, "AGE"),       
  typeCategorical = c("GENDER", "CENTER_C") 
)
metadeconfound_posthoc_adj <- MetaDeconfound(
  featureMat = drug_metabolite,
  metaMat = drug_metadata,
  QValues = metadeconfound_naive_adj$Qs,
  DValues = metadeconfound_naive_adj$Ds,
  typeContinuous = c(response, "AGE"),
  typeCategorical = c("GENDER", "CENTER_C")
)
adj_metadeconfound_df <- as.data.frame(cbind(metadeconfound_naive_adj$Ps,metadeconfound_naive_adj$Qs,metadeconfound_naive_adj$Ds,metadeconfound_posthoc_adj$status))
write.csv(adj_metadeconfound_df, "adj_metadeconfound.csv", row.names = TRUE)
save(adj_metadeconfound_df, file= "adj_metadeconfound.Rdata")

