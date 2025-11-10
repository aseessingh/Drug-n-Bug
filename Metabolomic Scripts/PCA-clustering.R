library(readr)
library(dplyr)
library(ropls)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
dose_serum_ids <- as.data.frame(read_csv("metadata3.csv",show_col_types = FALSE))
group_3 <- as.data.frame(read_csv("group 3.csv",show_col_types = FALSE))
AG <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))
metabolites <- as.data.frame(read_csv("MC_metabolites_wo_Drugs.csv",show_col_types = FALSE))
chemical_annotation <- as.data.frame(read_csv("chemical annotation.csv",show_col_types = FALSE))
#Preprocessing
AG <- AG %>% mutate(across(2:1359,log10))
metabolites <- metabolites %>% mutate(across(2:1359,log10)) %>% select(-`1,5-anhydroglucitol (1,5-AG)`)
group_3 <- group_3 %>% select(SampleID,GLYCATHB)
AG <- AG %>% select(MGS_PID,'1,5-anhydroglucitol (1,5-AG)') %>% rename(AG='1,5-anhydroglucitol (1,5-AG)')
dose_serum <- dose_serum_ids
dose_serum <- left_join(dose_serum,AG,by=c("SampleID"="MGS_PID"))
dose_serum <- left_join(dose_serum,group_3,by=c("SampleID"))
dose_serum <- na.omit(dose_serum)
dose_serum <- dose_serum %>% select(SampleID,GENDER,CENTER_C,everything())
pca <- opls(dose_serum[,7:8],predI = 2,orthoI = 0)
set.seed(123)
k_means <- kmeans(dose_serum[,7:8],centers = 2,nstart=100,iter.max = 1000)
dose_serum$Cluster <- as.factor(k_means$cluster)
dose_serum_plot <-  ggplot(dose_serum, aes(x = GLYCATHB, y = AG, color = Cluster)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  labs(title = "Scatter Plot with Clusters and 95% CI Ellipse",
       x = "GLYCATHB",
       y = "AG",
       color = "Cluster") +
  theme_minimal()
dose_serum_cluster_1 <- dose_serum %>% filter(Cluster == 1)
dose_serum_cluster_2 <- dose_serum %>% filter(Cluster == 2)
#Metabolite PLS-DA
dose_serum_cluster_1_metabolites <- metabolites %>% filter(MGS_PID %in% dose_serum_cluster_1$SampleID)
dose_serum_cluster_2_metabolites <- metabolites %>% filter(MGS_PID %in% dose_serum_cluster_2$SampleID)
dose_serum_cluster_1_metabolites$Cluster <- 1
dose_serum_cluster_2_metabolites$Cluster <- 2
dose_serum_metabolites <- rbind(dose_serum_cluster_1_metabolites,dose_serum_cluster_2_metabolites)
dose_serum_metabolites <- dose_serum_metabolites %>% select(MGS_PID,Cluster,everything())
dose_serum_pls <- opls(x=dose_serum_metabolites[,3:1359],y=dose_serum_metabolites$Cluster,predI = 3,orthoI = 0)
dose_serum_opls <- opls(x=dose_serum_metabolites[,3:1359],y=dose_serum_metabolites$Cluster,predI = 1,orthoI = 1)

#CA-PLS-DA
dose_serum_metabolites$Cluster <- NULL
combined_ca_pls_da <- left_join(dose_serum_metabolites,dose_serum,by=c("MGS_PID"="SampleID"))
extract_multiple_residuals <- function(data, metabolite_cols) {
  residuals_df <- data.frame(ID = 1:nrow(data))  # Optional: add ID for tracking
  
  for (col_num in metabolite_cols) {
    metabolite_name <- colnames(data)[col_num]
    # Wrap the metabolite name in backticks to handle spaces/special characters
    formula <- as.formula(paste0("`", metabolite_name, "` ~ AGE + GENDER + CENTER_C"))
    model <- lm(formula, data = data)
    residuals_df[[metabolite_name]] <- resid(model)
  }
  
  return(residuals_df)
}
ca_pls_x <- extract_multiple_residuals(data=combined_ca_pls_da,metabolite_cols = 2:1358)
ca_pls <- opls(x=ca_pls_x[,2:1358],y=dose_serum$Cluster,predI = 4,orthoI = 0)

#Loadings for All
pca_scores <- getScoreMN(pca)
pca_loadings <- getLoadingMN(pca)
pls_scores <- getScoreMN(dose_serum_pls)
pls_loadings <- getLoadingMN(dose_serum_pls)
opls_scores <- getScoreMN(dose_serum_opls)
opls_loadings <- getLoadingMN(dose_serum_opls)
opls_ortho_loadings <- getLoadingMN(dose_serum_opls,orthoL = TRUE)
ca_pls_scores <- getScoreMN(ca_pls)
ca_pls_loadings <- getLoadingMN(ca_pls)

#Function to Extract Metabolites from Loadings
combine_top_loadings <- function(mat1, mat2, mat3) {
  # Step 1: Convert matrices to data frames and preserve rownames
  df1 <- as.data.frame(mat1)
  df1$metabolite <- rownames(df1)
  
  df2 <- as.data.frame(mat2)
  df2$metabolite <- rownames(df2)
  
  df3 <- as.data.frame(mat3)  # OPLS orthogonal
  df3$metabolite <- rownames(df3)
  
  # Step 2: Get top 10 absolute values per column
  get_top_abs <- function(df, prefix, label, max_components = NULL) {
    n_components <- if (is.null(max_components)) ncol(df) - 1 else min(max_components, ncol(df) - 1)
    top_list <- lapply(seq_len(n_components), function(i) {
      col <- df[[i]]
      names(col) <- df$metabolite
      top_vals <- sort(abs(col), decreasing = TRUE)[1:10]
      top_names <- names(top_vals)
      data.frame(
        metabolite = top_names,
        value = col[top_names],
        component = paste0(prefix, "_", label, i),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, top_list)
  }
  
  top1 <- get_top_abs(df1, "PLS", "p")
  top2 <- get_top_abs(df2, "OPLS", "p")
  top3 <- get_top_abs(df3, "OPLS", "o")
  
  # Step 3: Combine all into one data frame
  combined <- rbind(top1, top2, top3)
  
  # Step 4: Reshape to wide format and merge duplicates
  library(tidyr)
  library(dplyr)
  
  wide_df <- combined %>%
    pivot_wider(names_from = component, values_from = value) %>%
    group_by(metabolite) %>%
    summarise(across(everything(), ~ if (all(is.na(.))) NA else first(na.omit(.))), .groups = "drop")
  
  return(wide_df)
}
top10df <- combine_top_loadings(pls_loadings,opls_loadings,opls_ortho_loadings)
chemical_annotation <- chemical_annotation %>% select(CHEMICAL_NAME,SUB_PATHWAY,SUPER_PATHWAY)
top10df <- left_join(top10df,chemical_annotation,by=c("metabolite"="CHEMICAL_NAME"))
top10df$SUB_PATHWAY[is.na(top10df$SUB_PATHWAY)] <- "Unknown"
top10df$SUPER_PATHWAY[is.na(top10df$SUPER_PATHWAY)] <- "Unknown"
top10df <- as.data.frame(top10df)
top10df[,2:6][is.na(top10df[,2:6])] <- 0
write.csv(top10df,"top 10 loadings with annotation.csv")

# Convert each variable to a data frame and write to CSV
write.csv(as.data.frame(pca_scores), file = "pca_scores.csv", row.names = TRUE)
write.csv(as.data.frame(pca_loadings), file = "pca_loadings.csv", row.names = TRUE)
write.csv(as.data.frame(pls_scores), file = "pls_scores.csv", row.names = TRUE)
write.csv(as.data.frame(pls_loadings), file = "pls_loadings.csv", row.names = TRUE)
write.csv(as.data.frame(opls_scores), file = "opls_scores.csv", row.names = TRUE)
write.csv(as.data.frame(opls_loadings), file = "opls_loadings.csv", row.names = TRUE)
write.csv(as.data.frame(opls_ortho_loadings), file = "opls_ortho_loadings.csv", row.names = TRUE)
write.csv(as.data.frame(ca_pls_scores), file = "ca_pls_scores.csv", row.names = TRUE)
write.csv(as.data.frame(ca_pls_loadings), file = "ca_pls_loadings.csv", row.names = TRUE)


#Export
write.csv(dose_serum," meta_clustering.csv")
