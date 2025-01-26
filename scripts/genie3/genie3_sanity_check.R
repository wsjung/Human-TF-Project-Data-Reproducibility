library(data.table)

# Need to update these path names. The oldest one is not actually the oldest
#old_path <- "/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_no_design/liver/genie3_post_processing/tf_by_gene.csv"
new_path <- "/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/liver/genie3_post_processing/tf_by_gene.csv"
#older_path <- "/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_08_08/liver/genie3_post_processing/tf_by_gene.csv"
oldest_path <-"/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/liver/old_genie3_post_processing/tf_by_gene.csv"

#old_file <- fread(old_path)
new_file <- fread(new_path)
#older_file <- fread(older_path)
oldest_file <- fread(oldest_path)

#old_file_good <- as.matrix(old_file, rownames=1)
new_file_good <- as.matrix(new_file, rownames=1)
#older_file_good <- as.matrix(older_file, rownames=1)
oldest_file_good <- as.matrix(oldest_file, rownames=1)

diff <- oldest_file_good - new_file_good
print(sum(abs(diff)))

#new_scores <- unlist(as.list(new_file_good))
#oldest_scores <- unlist(as.list(oldest_file_good))

#cor(oldest_scores, new_scores)

#cor.test(x=oldest_scores, y=new_scores, method = 'spearman')

#ranked_new_scores <- rank(-new_scores)
#ranked_oldest_scores <- rank(-oldest_scores)

#cor(ranked_oldest_scores, ranked_new_scores)

#plot(ranked_oldest_scores, ranked_new_scores, main='genie3 liver run ranks', xlab='run 1', ylab='run 2')

#diff2 = oldest_file_good - new_file_good
#print(sum(abs(diff2)))

#raw_ranks <- which(new_file_good >= tail(sort(new_file_good), n=10)[1], arr.ind=TRUE, useNames=TRUE)
#best_rows <- rownames(new_file_good)[raw_ranks[,1]]
#best_cols <- colnames(new_file_good)[raw_ranks[,2]]
#best_row_col_pairs <- rbind(best_rows,best_cols)

#new_file_good[best_row_col_pairs]
#ranks <- paste(rownames(new_file_good)[raw_ranks[,1]], colnames(new_file_good)[raw_ranks[,2]], new_file_good[raw_ranks[,1],raw_ranks[,2]], sep=',')
