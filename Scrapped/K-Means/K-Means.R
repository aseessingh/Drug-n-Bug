library(dplyr)
library(readr)
df <- as.data.frame(read_csv("combined_df_new.csv"))
df[,1] <- NULL
find_optimal_k <- function(df, cols, max_k = 10, scale_data = TRUE) {
  library(factoextra)
  
  # Subset and optionally scale the data
  data <- df[, cols]
  if (scale_data) {
    data <- scale(data)
  }
  
  # Elbow Method
  print("Elbow Method:")
  elbow_plot <- fviz_nbclust(data, kmeans, method = "wss", k.max = max_k) +
    ggtitle("Elbow Method for Optimal Clusters")
  print(elbow_plot)
  
  # Silhouette Method (factoextra uses its own silhouette implementation)
  print("Silhouette Method:")
  silhouette_plot <- fviz_nbclust(data, kmeans, method = "silhouette", k.max = max_k) +
    ggtitle("Silhouette Method for Optimal Clusters")
  print(silhouette_plot)
}

find_optimal_k(df,cols=12:254)
k_means <- kmeans(df[,12:254],centers = 2)
df$Cluster <- k_means$cluster
df <- df %>% select(SampleID,AGE,GENDERMale,GENDERFemale,CENTER_CDanemark,CENTER_CFrance,CENTER_CGermany,DOSAGE_METFORMIN_C,metformin,GLYCATHB,AG,Cluster,everything())
write_csv(df,"clustered_data.csv")
