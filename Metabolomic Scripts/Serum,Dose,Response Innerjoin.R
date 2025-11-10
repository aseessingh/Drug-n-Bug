library(dplyr)
library(readr)

serum_vs_metabolite <- read_csv("serum vs metabolite.csv", show_col_types = FALSE)
response_vs_metabolite <- read_csv("response vs metabolite.csv", show_col_types = FALSE)
dose_vs_metabolite <- read_csv("dose vs metab.csv",show_col_types = FALSE)
integrated_metabolite <- serum_vs_metabolite %>%
  inner_join(response_vs_metabolite, by=Name)%>%
  inner_join(dose_vs_metabolite,by=Name)

write.csv(integrated_metabolite, "integrated_metabolites.csv")