library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(DESeq2)
library(tidyverse)

#one read meta data 
metadata = read.csv("./meta.data.csv",header = TRUE, row.names = 1)

#read count matrix
counts <- read.csv("merged_counts.csv", header = TRUE)
rownames(counts) <- make.unique(as.character(counts[,1]))
counts <- counts[,-1]  # Remove the column now used as row names

base::colnames(counts) = row.names(metadata)



# Keep miRNAs with counts > 5 in at least 2 samples
keep <- rowSums(counts > 5) >= 2
counts_filtered <- counts[keep, ]





dds <- DESeqDataSetFromMatrix(countData = counts,  
                              colData = metadata,
                              design = ~ condition)  

dds <- DESeq(dds)
res <- results(dds)

view(res)

res_df <- as.data.frame(res)
res_filtered <- res_df[res_df$padj < 0.05 & !is.na(res_df$padj), ]






