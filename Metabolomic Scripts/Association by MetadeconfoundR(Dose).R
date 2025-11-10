library(dplyr)
library(metadeconfoundR)
library(readr)

#Input
metabolite_data_raw <- read_csv("MC_raw_metabolites.csv", show_col_types = FALSE)
treated_group <- read_csv("group 3.csv", show_col_types = FALSE)
healthy_group <- read_csv("group 8.csv", show_col_types = FALSE)
dosage_data <- read_csv("MC_Drug_doses.csv",show_col_types = FALSE)
dose_data <- c("DOSAGE_METFORMIN_C")
drug_info <- c("METFORMIN_C")
metabolite_data_without_drugs <- read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE)
drug_metabolite_in_data <- c("metformin")
#Make Code Fit
group3 <- treated_group
group8 <- healthy_group
drug_info_metadata <- drug_info
drug<-dose_data
#Filtering for Drug and Cleaning
drug_group3 <- group3 %>% filter(drug_info_metadata != 0) %>% select(all_of(c("SampleID","AGE", "GENDER", "CENTER_C")))
dosage_drug <- dosage_data %>% select(all_of(c("MCid", drug)))
drug_metadata <- left_join(drug_group3, dosage_drug, by = c("SampleID" = "MCid"))
drug_metadata <- drug_metadata %>% filter(rowSums(is.na(.) | . == "NA") == 0) %>% distinct(SampleID, .keep_all = TRUE)

#Extract matching data from metabolite dataset
drug_metabolite <- metabolite_data_without_drugs %>% filter(MGS_PID %in% drug_metadata[[1]])
drug_metabolite <- drug_metabolite %>% distinct(MGS_PID, .keep_all = TRUE)


#Add Group 8 to datasets for control
group8_meta <- group8 %>% select((all_of(c("SampleID","AGE", "GENDER", "CENTER_C"))))
group8_metabolites <- metabolite_data_without_drugs %>% filter(MGS_PID %in% group8_meta$SampleID)
group8_meta <- left_join(group8_meta,dosage_drug,by=c("SampleID"="MCid"))
drug_metadata$Status <- 1
group8_meta$Status <- 0
drug_metadata <- rbind(group8_meta, drug_metadata)
drug_metabolite <- rbind(drug_metabolite, group8_metabolites)

#Giving Row Names and harmonizing number in both datasets
drug_metadata <- as.data.frame(drug_metadata)
drug_metabolite <- as.data.frame(drug_metabolite)
drug_metadata <- drug_metadata %>% select(SampleID, everything())
drug_metadata <- drug_metadata %>% distinct(drug_metadata$SampleID, .keep_all=TRUE)
drug_metabolite <- drug_metabolite %>% distinct(drug_metabolite$MGS_PID, .keep_all = TRUE)
drug_metadata$`drug_metadata$SampleID` <- NULL
drug_metabolite$`drug_metabolite$MGS_PID` <- NULL
drug_metadata <- drug_metadata %>% filter(SampleID %in% drug_metabolite$MGS_PID)
rownames(drug_metadata) <- drug_metadata[, 1]
rownames(drug_metabolite) <- drug_metabolite[, 1]
drug_metadata[, 1] <- NULL
drug_metabolite[, 1] <- NULL
drug_metabolite <- drug_metabolite[match(rownames(drug_metadata),rownames(drug_metabolite)),,drop=FALSE]
write.csv(drug_metadata, "drug_metadata.csv",row.names=TRUE)
write.csv(drug_metabolite, "drug_metabolite.csv",row.names=TRUE)

#Bioavailability Delta Calculation
drug_consumers <- left_join(drug_group3, dosage_drug, by = c("SampleID" = "MCid"))
serum_drug_data <- metabolite_data_raw %>% select("MGS_PID",drug_metabolite_in_data)
drug_dose_serum <- left_join(drug_consumers,serum_drug_data,by=c("SampleID"="MGS_PID"))

#Metadeconfound Unadjusted for Demographics
drug_metadata_undaj <- drug_metadata %>% select(Status,drug)
metadeconfound_naive <- MetaDeconfound(featureMat = drug_metabolite, metaMat = drug_metadata_undaj, adjustMethod = "fdr", startStop = "naiveStop", typeContinuous = drug)
metadeconfound_posthoc <- MetaDeconfound(featureMat = drug_metabolite, metaMat = drug_metadata_undaj, QValues = metadeconfound_naive$Qs, DValues = metadeconfound_naive$Ds, typeContinuous = drug)
unadj_metadeconfound_df <- as.data.frame(cbind(metadeconfound_naive$Ps,metadeconfound_naive$Qs,metadeconfound_naive$Ds,metadeconfound_posthoc$status))
write.csv(unadj_metadeconfound_df, "unadj_metadeconfound.csv", row.names=TRUE)
save(unadj_metadeconfound_df, file= "unadj_metadeconfound.Rdata")

#Metadeconfound Adjusted for Demographics
drug_metadata <- drug_metadata %>% select(Status,everything())
metadeconfound_naive_adj <- MetaDeconfound(
  featureMat = drug_metabolite,
  metaMat = drug_metadata,
  adjustMethod = "fdr",
  startStop = "naiveStop",
  typeContinuous = c(drug, "AGE"),       
  typeCategorical = c("GENDER", "CENTER_C") 
)
metadeconfound_posthoc_adj <- MetaDeconfound(
  featureMat = drug_metabolite,
  metaMat = drug_metadata,
  QValues = metadeconfound_naive_adj$Qs,
  DValues = metadeconfound_naive_adj$Ds,
  typeContinuous = c(drug, "AGE"),
  typeCategorical = c("GENDER", "CENTER_C")
)
adj_metadeconfound_df <- as.data.frame(cbind(metadeconfound_naive_adj$Ps,metadeconfound_naive_adj$Qs,metadeconfound_naive_adj$Ds,metadeconfound_posthoc_adj$status))
write.csv(adj_metadeconfound_df, "adj_metadeconfound.csv", row.names = TRUE)
save(adj_metadeconfound_df, file= "adj_metadeconfound.Rdata")
write_csv(drug_dose_serum, "bioavailability_calc.csv")
