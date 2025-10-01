library('edgeR')
library('DESeq2')
library(magrittr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 4){ stop("Need exactly 4 args") }

tissues_tpm_dir <- args[1]
tissues_gene_counts_dir <- args[2]
tf_ids_path <- args[3]
output_dir <- args[4]

tf_ids <- read.csv(tf_ids_path, sep=',', header=TRUE)
names(tf_ids)<-str_replace_all(names(tf_ids), c(" " = "."))
tf_ensg_ids <- tf_ids['Ensembl.ID']

tissue_tpm_dirs <- list.dirs(tissues_tpm_dir, recursive=FALSE) # Obtain the paths to the inputs for each tissue
i <- 0
tpm_num_cols <- 0
for(tissue_tpm_dir in tissue_tpm_dirs) { 
	tissue <- basename(tissue_tpm_dir) # Grab tissue name
	print(tissue)
	if(i == 0){
		tpm_concat <- read.csv(file.path(tissue_tpm_dir, 'tpm.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)	
		tpm_num_cols <- tpm_num_cols + ncol(tpm_concat)
	}
	else{
		tpm <- read.csv(file.path(tissue_tpm_dir, 'tpm.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)
		tpm_num_cols <- tpm_num_cols + ncol(tpm)
		tpm_concat <- cbind(tpm_concat, tpm)
	}
	i <- i + 1
}
print(paste0('Number of tissues concatenated (tpm): ', i))
print(paste0('Expected number of columns in final tpm df: ', tpm_num_cols))
print(paste0('Number of columns in final tpm df: ', ncol(tpm_concat)))

tissue_gene_counts_dirs <- list.dirs(tissues_gene_counts_dir, recursive=FALSE) # Obtain the paths to the inputs for each tissue
j <- 0
gene_counts_num_cols <- 0
for(tissue_gene_counts_dir in tissue_gene_counts_dirs) {
        tissue <- basename(tissue_gene_counts_dir) # Grab tissue name
        print(tissue)
	if(j == 0){
                gene_counts_concat <- read.csv(file.path(tissue_gene_counts_dir, 'gene_counts.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)
                gene_counts_num_cols <- gene_counts_num_cols + ncol(gene_counts_concat)
        }
        else{
                gene_counts <- read.csv(file.path(tissue_gene_counts_dir, 'gene_counts.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)
                gene_counts_num_cols <- gene_counts_num_cols + ncol(gene_counts)
                gene_counts_concat <- cbind(gene_counts_concat, gene_counts)
        }
        j <- j + 1
}
print(paste0('Number of tissues concatenated (gene_counts): ', j))
print(paste0('Expected number of columns in final gene_counts df: ', gene_counts_num_cols))
print(paste0('Number of columns in final gene_counts df: ', ncol(gene_counts_concat)))

# Filter out low expression genes, such that CPM (counts per million) must be > 3 in >= 1.5% of a gene's samples
samplesThreshold <- floor(ncol(gene_counts_concat)*.015)
geneCountsCPM <- cpm(as.matrix(gene_counts_concat))
rowSumsInput <- geneCountsCPM > 3
rowSumsOutput <- rowSums(rowSumsInput)
expression_filter <- rowSumsOutput >= samplesThreshold
tpm_concat_filtered <- tpm_concat[expression_filter,]
print("gene_counts_concat")
print(head(gene_counts_concat, n = c(5, 5)))
print("samplesThreshold")
print(samplesThreshold)
print("geneCountsCPM")
print(head(geneCountsCPM, n=c(5,5)))
print("rowSumsInput")
print(head(rowSumsInput, n=c(5,5)))
print("rowSumsOutput")
print(head(rowSumsOutput))
print("expression_filter")
print(head(expression_filter))
print("tpm_concat")
print(head(tpm_concat, n = c(5, 5)))
print("tpm_concat_filtered")
print(head(tpm_concat_filtered, n=c(5,5)))
write.table(tpm_concat_filtered, file.path(output_dir, "np3_tpm_concat.txt"), sep='\t', row.names=TRUE, col.names=NA, quote=FALSE)

tf_ensg_ids <- unlist(tf_ensg_ids, use.names=FALSE)
tpm_genes <- row.names(tpm_concat_filtered)
head(tpm_genes)
# Get stable id for each ensg id in tpm_genes
tpm_genes_stable <- str_split_i(tpm_genes, "[.]", 1)
head(tpm_genes_stable)
tfs_filter <- tpm_genes_stable %in% tf_ensg_ids # Create a filter showing which tpm genes are tf's
tpm_tfs <- tpm_genes[tfs_filter]
print(paste0("Number of tfs: ", length(tpm_tfs)))
tf_tpm_concat_filtered <- tpm_concat_filtered[tfs_filter,]
print(paste0("Number of rows in tf_tpm_concat_filtered (should be # of tfs): ", nrow(tf_tpm_concat_filtered)))
write.table(tf_tpm_concat_filtered, file.path(output_dir, 'np3_tpm_concat_tfs.txt'), sep='\t', row.names=TRUE, col.names=NA, quote=FALSE)

print("done")

