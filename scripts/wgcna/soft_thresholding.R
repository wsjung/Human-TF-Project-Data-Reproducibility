.libPaths(c(.libPaths(),"/usr/local/lib/R/site-library"))
# Load WGCNA
library(WGCNA)

# Command-line arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]
beta <- as.numeric(args[3])

# Input data - a csv file with genes as columns and no other columns (e.g. sample column)
data <- read.csv(input_file, header = TRUE)

# Calculate the adjacency matrix using the selected power and pairwise complete observations - beta values are defined in the shell script
adjacency_matrix <- adjacency(data, power = beta, type = "unsigned", corOptions = list(use = "pairwise.complete.obs", method = "pearson"))

# Connectivity for each gene
connectivity <- colSums(adjacency_matrix) - 1  # Subtract 1 to exclude self-connection

# Mean connectivity
mean_connectivity <- mean(connectivity)

# Scale-free topology fit index
fit <- scaleFreeFitIndex(connectivity)

# Data frame with the results of each iteration - goes to intermediate output path
results <- data.frame(
  beta = beta,
  scale_free_fit_index = fit$Rsquared,
  mean_connectivity = mean_connectivity
)

# Save the df as a CSV file in the intermediate output path
write.csv(results, output_file, row.names = FALSE)

cat("Soft-thresholding iteration results saved to", output_file)
