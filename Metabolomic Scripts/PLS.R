library(readr)
library(dplyr)
library(ropls)
library(ggplot2)
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
set.seed(123)
k_means <- kmeans(dose_serum[,7:8],centers = 2,nstart=100)
dose_serum$Cluster <- as.factor(k_means$cluster)
dose_serum_cluster_1 <- dose_serum %>% filter(Cluster == 1)
dose_serum_cluster_2 <- dose_serum %>% filter(Cluster == 2)
dose_serum_cluster_1_metabolites <- metabolites %>% filter(MGS_PID %in% dose_serum_cluster_1$SampleID)
dose_serum_cluster_2_metabolites <- metabolites %>% filter(MGS_PID %in% dose_serum_cluster_2$SampleID)
dose_serum_cluster_1_metabolites$Cluster <- 1
dose_serum_cluster_2_metabolites$Cluster <- 2
dose_serum_metabolites <- rbind(dose_serum_cluster_1_metabolites,dose_serum_cluster_2_metabolites)
dose_serum_metabolites <- dose_serum_metabolites %>% select(MGS_PID,Cluster,everything())

#PLS-DA
pls <- opls(x=dose_serum_metabolites[,2:1358],y=dose_serum$Cluster,predI=3,orthoI=0,permI = 1000)
opls <- opls(x=dose_serum_metabolites[,2:1358],y=dose_serum$Cluster,predI=1,orthoI=1,permI = 1000)
#W,W' and L
pls_w <- as.data.frame(getWeightMN(pls))
pls_w_prime <- as.data.frame(pls@weightStarMN)
pls_l <- as.data.frame(getLoadingMN(pls))
pls_scores <- as.data.frame(getScoreMN(pls))
opls_w <- as.data.frame(getWeightMN(opls))
opls_w_ortho <- as.data.frame(getWeightMN(opls, orthoL = TRUE))
opls_w_prime <- as.data.frame(opls@weightStarMN)
opls_l <- as.data.frame(getLoadingMN(opls))
opls_l_ortho <- as.data.frame(getLoadingMN(opls, orthoL = TRUE))
opls_scores <- as.data.frame(getScoreMN(opls))
opls_scores_ortho <- as.data.frame(getScoreMN(opls,orthoL=TRUE))
opls_scores <- cbind(opls_scores,opls_scores_ortho)
dose_serum_cluster <- as.data.frame(dose_serum$Cluster)
opls_scores <- cbind(opls_scores,dose_serum_cluster)
pls_scores <- cbind(pls_scores,dose_serum_cluster)
#Extract Top 10
combine_top_loadings_general <- function(named_dataframes, top_n = 10, max_components = NULL) {
  library(dplyr)
  library(tidyr)
  
  get_top_abs <- function(df, label, max_components = NULL) {
    df$metabolite <- rownames(df)
    n_components <- if (is.null(max_components)) ncol(df) - 1 else min(max_components, ncol(df) - 1)
    
    top_list <- lapply(seq_len(n_components), function(i) {
      col <- df[[i]]
      names(col) <- df$metabolite
      top_vals <- sort(abs(col), decreasing = TRUE)[1:top_n]
      top_names <- names(top_vals)
      data.frame(
        metabolite = top_names,
        value = col[top_names],
        component = paste0(label, "_", i),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, top_list)
  }
  
  # Apply to all data frames
  all_top <- do.call(rbind, lapply(names(named_dataframes), function(name) {
    get_top_abs(named_dataframes[[name]], name, max_components)
  }))
  
  # Reshape to wide format
  wide_df <- all_top %>%
    pivot_wider(names_from = component, values_from = value) %>%
    group_by(metabolite) %>%
    summarise(across(everything(), ~ if (all(is.na(.))) NA else first(na.omit(.))), .groups = "drop")
  
  return(wide_df)
}
named_dataframes <- list(
  PLS_w = pls_w,
  PLS_w_prime = pls_w_prime,
  PLS_l = pls_l,
  OPLS_w = opls_w,
  OPLS_w_ortho = opls_w_ortho,
  OPLS_w_prime = opls_w_prime,
  OPLS_l = opls_l,
  OPLS_l_ortho = opls_l_ortho
)
result_df <- combine_top_loadings_general(named_dataframes)
chemical_annotation <- chemical_annotation %>% select(CHEMICAL_NAME,SUB_PATHWAY,SUPER_PATHWAY)
result_df <- left_join(result_df,chemical_annotation, by=c("metabolite"="CHEMICAL_NAME"))

#t-test along p1
run_group_ttest_by_index <- function(df, group_col_index, value_col_index) {
  # Check if indices are valid
  if (group_col_index > ncol(df) || value_col_index > ncol(df)) {
    stop("One or both column indices are out of bounds.")
  }
  
  # Extract column names from indices
  group_col <- names(df)[group_col_index]
  value_col <- names(df)[value_col_index]
  
  # Ensure the grouping column is a factor
  df[[group_col]] <- as.factor(df[[group_col]])
  
  # Check that there are exactly two levels
  if (length(levels(df[[group_col]])) != 2) {
    stop("Grouping column must have exactly two levels for a t-test.")
  }
  
  # Perform the t-test
  formula <- as.formula(paste(value_col, "~", group_col))
  t_test_result <- t.test(formula, data = df)
  
  return(t_test_result)
}
export_ttest_results_to_csv <- function(ttest_result, filename = "ttest_results.csv") {
  # Create a data frame with relevant statistics
  result_df <- data.frame(
    statistic = ttest_result$statistic,
    p_value = ttest_result$p.value,
    conf_int_low = ttest_result$conf.int[1],
    conf_int_high = ttest_result$conf.int[2],
    mean_group1 = ttest_result$estimate[1],
    mean_group2 = ttest_result$estimate[2]
  )
  
  # Write to CSV
  write.csv(result_df, file = filename, row.names = FALSE)
}

pls_t <- run_group_ttest_by_index(df=pls_scores,group_col = 4,value_col = 2)
opls_t <- run_group_ttest_by_index(df=opls_scores,group_col = 3,value_col = 1)

#Correlation

#Export
write.csv(result_df,"all parameters of pls.csv")
write.csv(pls_scores,"pls scores.csv")
write.csv(opls_scores, "opls scores.csv")
export_ttest_results_to_csv(pls_t,"pls_t.csv")
export_ttest_results_to_csv(opls_t, "opls_t.csv")
save(pls, file = "pls.Rdata")
save(opls, file = "opls.Rdata")
opls_w <- cbind(opls_w,opls_w_ortho)
write.csv(opls_w,"opls weights.csv")
write.csv(pls_w,"pls weight.csv")
write.csv(pls_w_prime,"pls weight'.csv")
write.csv(opls_w_prime,"opls weight'.csv")
write.csv(pls_l,"pls loadings.csv")
opls_l <- cbind(opls_l,opls_l_ortho)
write.csv(opls_l,"opls loadings.csv")
