# Load necessary libraries
library(WGCNA)

# Command-line arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
adjacency_output_file <- args[2]
tom_output_file <- args[3]
diss_tom_output_file <- args[4]
beta <- as.numeric(args[5])

# Input data
data <- read.csv(input_file, header = TRUE)

# Extract gene IDs from column headers for later use when constructing the matrixes
gene_ids <- colnames(data)

# Adjacency matrix using the selected power (best beta value is obtained from the soft_thresholding script)
adjacency_matrix <- adjacency(data, power = beta)
rownames(adjacency_matrix) <- gene_ids
colnames(adjacency_matrix) <- gene_ids

# Verification that the adjacency matrix has correct gene IDs
cat("Adjacency matrix first 5 gene IDs (row):", head(rownames(adjacency_matrix)), "\n")
cat("Adjacency matrix first 5 gene IDs (column):", head(colnames(adjacency_matrix)), "\n")

# Connectivity for each gene
connectivity <- colSums(adjacency_matrix) - 1  # Subtract 1 to exclude self-connection

# Mean connectivity
mean_connectivity <- mean(connectivity)

# Scale-free topology fit index
fit <- scaleFreeFitIndex(connectivity)

# Print the scale-free topology fit index and mean connectivity (should be same as the result of this beta from the first script)
cat("Scale-Free Topology Fit Index for beta =", beta, ":", fit$Rsquared, "\n")
cat("Mean Connectivity for beta =", beta, ":", mean_connectivity, "\n")

# Topological Overlap Matrix (TOM)
tom <- TOMsimilarity(adjacency_matrix)

# Verify the dimensions of TOM matrix
if (all(dim(tom) == dim(adjacency_matrix))) {
  cat("TOM matrix dimensions match the adjacency matrix dimensions.\n")
} else {
  stop("Mismatch in TOM matrix dimensions.\n")
}

# Ensuring the gene IDs are correctly assigned to the TOM matrix
rownames(tom) <- gene_ids
colnames(tom) <- gene_ids

# Verify the assignment by checking if the order of gene IDs matches the original
if (all(rownames(tom) == gene_ids) && all(colnames(tom) == gene_ids)) {
  cat("Gene IDs are correctly assigned to the TOM matrix.\n")
} else {
  stop("Mismatch in gene IDs assignment to the TOM matrix.\n")
}

# The dissimilarity TOM (1 - TOM)
diss_tom <- 1 - tom

# Ensure the gene IDs are correctly assigned to the dissimilarity TOM matrix
rownames(diss_tom) <- gene_ids
colnames(diss_tom) <- gene_ids

# Verify the assignment by checking if the order of gene IDs matches the original
if (all(rownames(diss_tom) == gene_ids) && all(colnames(diss_tom) == gene_ids)) {
  cat("Gene IDs are correctly assigned to the dissimilarity TOM matrix.\n")
} else {
  stop("Mismatch in gene IDs assignment to the dissimilarity TOM matrix.\n")
}

# Save the TOM as a CSV file
write.csv(tom, tom_output_file, row.names = TRUE, quote = FALSE)

# Save the dissimilarity TOM as a CSV file
write.csv(diss_tom, diss_tom_output_file, row.names = TRUE, quote = FALSE)

cat("TOM saved to", tom_output_file, "\n")
cat("Dissimilarity TOM saved to", diss_tom_output_file, "\n")
