library(readr)

unadj <- read_csv("",show_col_types = FALSE)
adj <- read.csv("",show_col_types = FALSE)

intersection <- intersect(length(unadj$Name),length(adj$Name))
union <- length(unadj$Name) + length(adj$Name) - intersection

jaccard <- intersection / union

print(jaccard)