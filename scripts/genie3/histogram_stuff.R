library(data.table)

genie3_path <- '/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07' # GENIE3 path for tissue

tissue <- 'liver'

print(tissue)

tissue_tf_gene_genie3 <- fread(file.path(genie3_path, tissue, 'genie3_post_processing', 'tf_by_gene.csv'))
tissue_tf_gene_genie3 <- as.matrix(tissue_tf_gene_genie3, rownames=1)

# Make sure there's no autoregulation (elements with same row and column name should be 0)
old_tissue_tf_gene_genie3 <- tissue_tf_gene_genie3
old_tissue_tf_gene_genie3[outer(rownames(old_tissue_tf_gene_genie3), colnames(old_tissue_tf_gene_genie3), "==")] <- 0
if(!identical(tissue_tf_gene_genie3, old_tissue_tf_gene_genie3)){
	stop('original matrix was not autoregulated')
}

# Ensure that values are from 0 to 1
print(paste0("min of tissue_tf_gene_genie3: ", min(tissue_tf_gene_genie3)))
print(paste0("max of tissue_tf_gene_genie3: ", max(tissue_tf_gene_genie3)))

# Create histogram and kde distribution of scores
hist(tissue_tf_gene_genie3, breaks=100, freq = FALSE, main = paste0(tissue, " TF x Gene GENIE3 Scores"))
dens <- density(tissue_tf_gene_genie3)
lines(dens)

