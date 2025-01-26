# Script for picking a threshold for WGCNA beta coefficient - much of this was taken from https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)

# Load dataset - should be vst transformed, genes as columns, no row names
tissue <- 'liver'

vst_output_path <- file.path('/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output', tissue, 'no_design_gene_counts_vst_wgcna_genie3.csv')
VST <- read.csv(vst_output_path, header=TRUE)

# Choose a set of soft-thresholding powers
powers = c(1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(
  VST,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
#  RsquaredCut = 0.85
  )

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste0(tissue, " Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste0(tissue, " Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

