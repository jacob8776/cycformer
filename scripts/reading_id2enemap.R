library(tidyverse)
library(ggplot2)

setwd("~/Desktop/jacob/mcycdb/MCycDB/")

metadata <- read_table("id2genemap.txt", col_names = F) 

colnames(metadata) <- c("id", "gene", "database")


metadata <- filter(metadata, is.na(database) == FALSE | is.na(id) == FALSE | is.na(gene) == FALSE)

write.csv(metadata, "metadata_filtered.csv", row.names = F)


data <- read_table("MCycDB_2021.faa", col_names = F)













