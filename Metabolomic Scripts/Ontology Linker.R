library(dplyr)
library(readr)

annotated <- read_csv("Chemical Annotations Integrated with logP.csv",show_col_types = FALSE)
ontology <- read_csv("Metabolite_Ontology_InHouse_08042025.csv",show_col_types = FALSE)

ontology <- ontology %>% select (CHEMICAL_NAME, ontology)
annotated_final <- inner_join(annotated, ontology, by=c("CHEMICAL_NAME"))
write_csv(annotated_final, "Chemical Annotations Integrated Final.csv")
