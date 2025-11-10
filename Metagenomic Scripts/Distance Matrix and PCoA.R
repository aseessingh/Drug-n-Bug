library(vegan)
library(readr)
library(dplyr)
library(ggplot2)
#Load preprocessed data
HC <- as.data.frame(read_csv("HC.csv",show_col_types = FALSE))
UMMC <- as.data.frame(read_csv("UMMC.csv", show_col_types = FALSE))
TREATMENT_GROUP <- as.data.frame(read_csv("TG_SU.csv",show_col_types = FALSE))

TG_MET <- TREATMENT_GROUP

#Change IDs to column names and export order for metadata df
rownames(HC) <- HC[,1]
rownames(UMMC) <- UMMC[,1]
rownames(TG_MET) <- TG_MET[,1]
vector_list <- list(
  HC = HC[,1],
  UMMC = UMMC[,1],
  Treatment = TG_MET[,1]
)
row_order <- bind_rows(
  lapply(names(vector_list), function(name) {
    data.frame(Column = name, IDs = vector_list[[name]])
  })
)
HC[,1] <- NULL
UMMC[,1] <- NULL
TG_MET[,1] <- NULL

#Distance Matrix and PcoA
HC_dist <- vegdist(HC, method = "bray")
UMMC_dist <- vegdist(UMMC, method = "bray")
TG_MET_dist <- vegdist(TG_MET, method="bray")
HC_PcoA <- cmdscale(HC_dist, k=2, eig=TRUE)
UMMC_PcoA <- cmdscale(UMMC_dist,k=2,eig=TRUE)
TG_MET_PcoA <- cmdscale(TG_MET_dist, k=2, eig=TRUE)

#Visualisation of PCoA
HC_df <- data.frame(
  pc1 = HC_PcoA$points[,1],
  pc2 = HC_PcoA$points[,2],
  Group = "Healthy"
)
UMMC_df <- data.frame(
  pc1 = UMMC_PcoA$points[,1],
  pc2= UMMC_PcoA$points[,2],
  Group = "Unhealthy"
)
TG_MET_df <- data.frame(
  pc1 = TG_MET_PcoA$points[,1],
  pc2 = TG_MET_PcoA$points[,2],
  Group = "Drug Treated"
)
PCoA_df <- rbind(HC_df,UMMC_df,TG_MET_df)
plot <- ggplot(PCoA_df, aes(x = pc1, y = pc2, colour = Group)) +
  geom_point(size = 2) +
  stat_ellipse(level= 0.95)+
  labs(
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    title = "Trends in Microbial Communities"
  ) +
  theme_minimal()

save(HC_dist, file= "HC_dist.Rdata")
save(UMMC_dist, file= "UMMC_dist.Rdata")
save(TG_MET_dist, file= "Gliclazide_dist.Rdata")
ggsave("glic.png", plot = plot, width = 10, height = 8, dpi = 600)

