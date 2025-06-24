
library(tidyverse)
library(dplyr)


# Set working directory
setwd("../Desktop/microRNA/")

# List CSV files
files <- list.files(path = "./expression matrix/", pattern = "*.csv", full.names = TRUE)


# Read the first file
df1 <- read.csv(files[1], sep = "\t", header = TRUE)
head(df1)
df1 <- df1[, c("X.miRNA", "read_count")]
colnames(df1) <- c("miRNA", tools::file_path_sans_ext(files[1]))

# Read and merge the rest
for (i in 2:length(files)) {
  df <- read.csv(files[i], sep = "\t", header = TRUE)
  df <- df[, c("X.miRNA", "read_count")]
  colnames(df) <- c("miRNA", tools::file_path_sans_ext(files[i]))
  df1 <- merge(df1, df, by = "miRNA", all = TRUE)
}

# Replace NA with 0
df1[is.na(df1)] <- 0

# Save merged count matrix
write.csv(df1, "merged_counts.csv", row.names = FALSE)












