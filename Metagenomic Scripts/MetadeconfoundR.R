library(readr)
library(tidyverse)
library(metadeconfoundR)

#Input Parameters
drug_intake <- as.data.frame(read_csv("metadata3.csv"))
non_drug_intake <- as.data.frame(read_csv("metadata4.csv"))
metagenome <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.mgsamples.v6.r.csv"))
drug_dose <- as.data.frame(read_csv("MC_Drug_doses.csv"))
drug <- c("metformin")
metagenome_annotation <- as.data.frame(read_csv("hub.taxon_adjusted.MGS.down.10000000.variables.v6.r.csv"))
#Preprocessing
metagenome$MGSampleID <- gsub("M0_","",metagenome$MGSampleID)
metagenome_annotation <- metagenome_annotation %>% dplyr::select(VariableID,DisplayName,Notes)
log_transform_normalize <- function(df) {
  # Exclude the first column
  data_to_transform <- df[, -1]
  
  # Apply log transformation (adding 1 to avoid log(0))
  log_transformed <- log(data_to_transform + 1)
  
  # Min-max normalization function
  normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Apply normalization
  normalized <- as.data.frame(lapply(log_transformed, normalize))
  
  # Combine the first column back with the transformed data
  result <- cbind(df[, 1, drop = FALSE], normalized)
  
  return(result)
}
metagenome <- log_transform_normalize(metagenome)
non_drug_intake <- non_drug_intake %>% filter(Status == 0)
drug_intake <- drug_intake %>% select(SampleID,AGE,GENDER,CENTER_C)
non_drug_intake <- dplyr::left_join(non_drug_intake,drug_dose,by=c("SampleID"="MCid"))
filter_non_drug_intake <- function(df, drug) {
  if (drug == "metformin") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_DPPIV_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_SU_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else if (drug == "su") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_DPPIV_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_METFORMIN_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else if (drug == "dppiv") {
    non_drug_intake_0 <- df %>%
      filter(DOSAGE_ACARBOSE_C == 0 & DOSAGE_BOLUS_C == 0 & DOSAGE_GLP_1_C == 0 &
               DOSAGE_INSULIN_C == 0 & DOSAGE_SGLT2_C == 0 & DOSAGE_SU_C == 0 & DOSAGE_METFORMIN_C == 0) %>%
      select(SampleID, AGE, GENDER, CENTER_C)
  } else {
    stop("Invalid drug name")
  }
  
  return(non_drug_intake_0)
}
non_drug_intake_0 <- filter_non_drug_intake(df=non_drug_intake,drug=drug)
non_drug_intake <- non_drug_intake %>% select(SampleID,AGE,GENDER,CENTER_C)
metagenome_drug <- metagenome %>% filter(MGSampleID %in% drug_intake$SampleID)
metagenome_non_drug <- metagenome %>% filter(MGSampleID %in% non_drug_intake$SampleID)
metagenome_non_drug_0 <- metagenome %>% filter (MGSampleID %in% non_drug_intake_0$SampleID)
rownames_function <- function(df) {
  rownames(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}
metagenome_drug <- rownames_function(metagenome_drug)
metagenome_non_drug <- rownames_function(metagenome_non_drug)
metagenome_non_drug_0 <- rownames_function(metagenome_non_drug_0)
drug_intake$Group <- 1
non_drug_intake$Group <- 0
non_drug_intake_0$Group <- 0
drug_intake <- drug_intake %>% select(SampleID,Group, everything())
non_drug_intake <- non_drug_intake %>% select(SampleID,Group, everything())
non_drug_intake_0 <- non_drug_intake_0 %>% select(SampleID,Group, everything())
drug_non_drug_md <- rbind(drug_intake,non_drug_intake)
drug_non_drug_0_md <- rbind(drug_intake,non_drug_intake_0)
drug_non_drug_mg <- rbind(metagenome_drug,metagenome_non_drug)
drug_non_drug_0_mg <- rbind (metagenome_drug,metagenome_non_drug_0)

#MetadeconfoundR
drug_non_drug_md <- drug_non_drug_md %>% distinct(SampleID,.keep_all = TRUE)
rownames(drug_non_drug_md) <- drug_non_drug_md[,1]
drug_non_drug_md[,1] <- NULL
drug_non_drug_md <- drug_non_drug_md[match(rownames(drug_non_drug_mg), rownames(drug_non_drug_md)), ]
drug_non_drug_metadeconfound_run<- MetaDeconfound(
  featureMat = drug_non_drug_mg,
  metaMat = drug_non_drug_md,
  typeCategorical = c("GENDER","CENTER_C"),
  typeContinuous = c("AGE"),
  QCutoff = 0.05,
  DCutoff = 0.2,
  nnodes = 16,
  adjustMethod = "fdr"
)
drug_non_drug_metadeconfound <- as.data.frame(cbind(drug_non_drug_metadeconfound_run$Ps,drug_non_drug_metadeconfound_run$Qs,drug_non_drug_metadeconfound_run$Ds,drug_non_drug_metadeconfound_run$status))
colnames(drug_non_drug_metadeconfound) <- make.unique(colnames(drug_non_drug_metadeconfound))
drug_non_drug_metadeconfound <- drug_non_drug_metadeconfound %>%
  filter(Group.3 %in% c("OK_nc", "OK_sd")) %>% select (Group,Group.1,Group.2) %>% rename(p = Group, q = Group.1 , Rho = Group.2)
drug_non_drug_metadeconfound$Variable <- rownames(drug_non_drug_metadeconfound)
drug_non_drug_metadeconfound <- left_join(drug_non_drug_metadeconfound,metagenome_annotation,by=c("Variable"="VariableID"))
rownames(drug_non_drug_metadeconfound) <- NULL
write.csv(drug_non_drug_metadeconfound,"drug_non_drug_metadeconfound.csv")
#
drug_non_drug_0_md <- drug_non_drug_0_md %>% distinct(SampleID,.keep_all = TRUE)
rownames(drug_non_drug_0_md) <- drug_non_drug_0_md[,1]
drug_non_drug_0_md[,1] <- NULL
drug_non_drug_0_md <- drug_non_drug_0_md[match(rownames(drug_non_drug_0_mg), rownames(drug_non_drug_0_md)), ]
drug_non_drug_0_metadeconfound_run<- MetaDeconfound(
  featureMat = drug_non_drug_0_mg,
  metaMat = drug_non_drug_0_md,
  typeCategorical = c("GENDER","CENTER_C"),
  typeContinuous = c("AGE"),
  QCutoff = 0.05,
  DCutoff = 0.2,
  nnodes = 16,
  adjustMethod = "fdr"
)
drug_non_drug_0_metadeconfound <- as.data.frame(cbind(drug_non_drug_0_metadeconfound_run$Ps,drug_non_drug_0_metadeconfound_run$Qs,drug_non_drug_0_metadeconfound_run$Ds,drug_non_drug_0_metadeconfound_run$status))
colnames(drug_non_drug_0_metadeconfound) <- make.unique(colnames(drug_non_drug_0_metadeconfound))
drug_non_drug_0_metadeconfound <- drug_non_drug_0_metadeconfound %>%
  filter(Group.3 %in% c("OK_nc", "OK_sd")) %>% select (Group,Group.1,Group.2) %>% rename(p = Group, q = Group.1 , Rho = Group.2)
drug_non_drug_0_metadeconfound$Variable <- rownames(drug_non_drug_0_metadeconfound)
drug_non_drug_0_metadeconfound <- left_join(drug_non_drug_0_metadeconfound,metagenome_annotation,by=c("Variable"="VariableID"))
rownames(drug_non_drug_0_metadeconfound) <- NULL
write.csv(drug_non_drug_0_metadeconfound,"drug_non_drug_0_metadeconfound.csv")

