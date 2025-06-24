BiocManager::install("edgeR")
library(edgeR)

#one read meta data 
metadata = read.csv("./meta.data.csv",header = TRUE, row.names = 1)

#read count matrix
counts <- read.csv("merged_counts.csv", header = TRUE)
rownames(counts) <- make.unique(as.character(counts[,1]))
counts <- counts[,-1]  # Remove the column now used as row names

base::colnames(counts) = row.names(metadata)



# Load counts and metadata
counts <- read.csv("merged_counts.csv", row.names = 1)
group <- factor(c("WT", "WT", "WT", "KO", "KO", "KO"))  # adjust based on your samples

# Create DGEList object
dge <- DGEList(counts = counts, group = group)

# Filter out lowly expressed miRNAs
keep <- filterByExpr(dge)
dge <- dge[keep,, keep.lib.sizes=FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Estimate dispersion
dge <- estimateDisp(dge)

# Fit model and test
fit <- exactTest(dge)

# Results table
res <- topTags(fit, n = Inf)
head(res$table)



