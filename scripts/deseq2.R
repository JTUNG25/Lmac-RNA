# scripts/deseq2.R

# Load the DESeq2 library
library(DESeq2)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
samples_file <- args[2]
contrast_str <- args[3]
output_file <- args[4]

# Read the counts data
countData <- read.table(counts_file, header=TRUE, row.names=1)

# Read the sample data
colData <- read.table(samples_file, header=TRUE, row.names=1)

# Make sure the column names of the count data match the row names of the sample data
# and are in the same order
countData <- countData[, rownames(colData)]

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Run the DESeq analysis
dds <- DESeq(dds)

# Parse the contrast string
contrast_parts <- unlist(strsplit(contrast_str, "_vs_"))
condition1 <- contrast_parts[1]
condition2 <- contrast_parts[2]

# Get the results for the specified contrast
res <- results(dds, contrast=c("condition", condition1, condition2))

# Write the results to a CSV file
write.csv(as.data.frame(res), file=output_file)
