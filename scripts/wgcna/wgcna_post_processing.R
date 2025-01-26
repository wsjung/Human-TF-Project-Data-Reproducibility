library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2){ stop("Need exactly 2 args") }

wgcna_output_path <- args[1]
tf_ids_path <- args[2]

tissue_TOM <- fread(file.path(wgcna_output_path, 'tom.csv'))
tissue_TOM <- as.matrix(tissue_TOM, rownames=1)                                                                                                                                                                                                                                                                         
# Ensure that autoregulation isn't allowed, by setting elements with same row and column name equal to 0
tissue_TOM[outer(rownames(tissue_TOM), colnames(tissue_TOM), "==")] <- 0                                                                                                                                                                                                                                                
tf_ids <- read.csv(tf_ids_path, sep='\t', header=TRUE)
tf_ensg_ids <- tf_ids['TF_ensg']                                                                                                                                                                                                                                                                                        
# Filter the tissue genes based on which are tf's.
tf_ensg_ids <- unlist(tf_ensg_ids, use.names=FALSE) # Make sure tissue genes and tf_ensg_ids are both character vectors, for easy comparison
tissue_genes <- row.names(tissue_TOM)
tissue_tfs_filter <- tissue_genes %in% tf_ensg_ids                                                                                                          
tissue_tfs <- tissue_genes[tissue_tfs_filter]
print(paste0("Number of tissue tfs: ", length(tissue_tfs)))                                                                                                                                                                                                                                                             
# Create TFs x genes correlation matrix for each tissue (MxN where M is # TFs, N is # genes)                                                                
tissue_tf_gene_TOM <- tissue_TOM[tissue_tfs_filter,]
print(paste0("Number of rows in tissue_tf_gene_TOM (should be # of tissue tfs): ", nrow(tissue_tf_gene_TOM)))
write.csv(tissue_tf_gene_TOM, file.path(wgcna_output_path, 'tf_gene_tom.csv'), quote=FALSE)                                                                                                                                                                                                                             
# Ensure that values are from 0 to 1                                                                                                                        
print(paste0("min of tissue_tf_gene_TOM: ", min(tissue_tf_gene_TOM)))
print(paste0("max of tissue_tf_gene_TOM: ", max(tissue_tf_gene_TOM)))

print('done') 


