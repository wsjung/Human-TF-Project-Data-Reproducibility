library('edgeR')
library('DESeq2')
library(magrittr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){ stop("Need exactly 3 args") }

tissue_tpm_dir <- args[1]
tissue_gene_counts_dir <- args[2]
tf_ids_path <- args[3]

tf_ids <- read.csv(tf_ids_path, sep=',', header=TRUE)
# replace spaces in column names with underlines
names(tf_ids)<-str_replace_all(names(tf_ids), c(" " = "."))
tf_ensg_ids <- tf_ids['Ensembl.ID']

tissue <- basename(tissue_tpm_dir) # Grab tissue name
print(tissue)
tpm <- read.csv(file.path(tissue_tpm_dir, 'tpm.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)
geneCounts <- read.csv(file.path(tissue_gene_counts_dir, 'gene_counts.txt'), sep='\t', row.names = 1, check.names=FALSE, header=TRUE)

# Filter out low expression genes, such that CPM (counts per million) must be > 3 in >= 1.5% of a gene's samples
samplesThreshold <- floor(ncol(geneCounts)*.015)
geneCountsCPM <- cpm(as.matrix(geneCounts))
rowSumsInput <- geneCountsCPM > 3
rowSumsOutput <- rowSums(rowSumsInput)
expression_filter <- rowSumsOutput >= samplesThreshold
tpmFiltered <- tpm[expression_filter,]
print("geneCounts")
print(head(geneCounts, n = c(5, 5)))
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
print("tpm")
print(head(tpm, n = c(5, 5)))
print("tpmFiltered")
print(head(tpmFiltered, n=c(5,5)))
write.table(tpmFiltered, file.path(tissue_tpm_dir, "np3_tpm.txt"), sep='\t', row.names=TRUE, col.names=NA, quote=FALSE)

# Create TFs x genes matrix for each tissue (MxN where M is # TFs, N is # genes)
tf_ensg_ids <- unlist(tf_ensg_ids, use.names=FALSE)
tpm_genes <- row.names(tpmFiltered)
head(tpm_genes)
# Get stable id for each ensg id in tpm_genes
tpm_genes_stable <- str_split_i(tpm_genes, "[.]", 1)
head(tpm_genes_stable)
tfs_filter <- tpm_genes_stable %in% tf_ensg_ids # Create a filter showing which tpm genes are tf's
tpm_tfs <- tpm_genes[tfs_filter]
print(paste0("Number of tfs: ", length(tpm_tfs)))
tf_tpmFiltered <- tpmFiltered[tfs_filter,]
print(paste0("Number of rows in tf_tpmFiltered (should be # of tfs): ", nrow(tf_tpmFiltered)))
write.table(tf_tpmFiltered, file.path(tissue_tpm_dir, 'np3_tpm_tfs.txt'), sep='\t', row.names=TRUE, col.names=NA, quote=FALSE)

print("done")

