# 02_gsea_annotation.R

# Load libraries
library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(tidyverse)

# Load DE results
res_df <- read.csv("deseq2_all_results.csv", row.names = 1)

# Remove NA values
res_df <- res_df[!is.na(res_df$padj), ]

# Create ranked gene list for GSEA (sorted by log2FoldChange)
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA using KEGG or GO — Arabidopsis-specific
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.At.tair.db,
                 ont = "BP",             # Biological Process
                 keyType = "TAIR",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

# Save GSEA results
write.csv(as.data.frame(gsea_go), file = "gsea_go_results.csv")

# Plot top pathways
dotplot(gsea_go, showCategory = 15) + ggtitle("GSEA - GO BP")
