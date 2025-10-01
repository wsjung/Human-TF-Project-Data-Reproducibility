library(glmnet)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
input_path <- args[1]
batch_num <- args[2]
output_path <- args[3]


# create output dir if not exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive=TRUE)
}

# load data (transpose expression for glmnet)
df_data <- read.table(file.path(input_path, "np3_tpm.txt"),
                      header=TRUE, sep="\t", row.names=1)
df_tfs <- read.table(file.path(input_path, "np3_tpm_tfs.txt"),
                     header=TRUE, sep="\t", row.names=1)
list_genes <- readLines(file.path(input_path, "bart_batches", paste0(batch_num,".txt")))
N_genes <- length(list_genes)

# preprocess data
df_data <- df_data %>% t
df_tfs <- df_tfs %>% t


# list to hold coefs
lasso_results <- list()

idx <- 1
for(gene in list_genes) {
  print(paste0(gene," [",idx,"/",N_genes,"]"))
  
  # subset gene
  vct_gene <- df_data[,gene]
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(df_tfs, vct_gene, alpha = 1)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  
  # find best model
  best_model <- glmnet(df_tfs, vct_gene, alpha = 1, lambda = best_lambda)
  
  # coefficients
  coefs <- as.vector(coef(best_model))
  names(coefs) <- rownames(coef(best_model))
  
  # store coefs
  lasso_results[[gene]] <- coefs
  
  # update index
  idx = idx + 1
}

# concatenate lasso results
df_lasso <- do.call(rbind, lasso_results)
df_lasso[is.na(df_lasso)] <- 0
df_lasso <- as.data.frame(df_lasso) %>% select(-`(Intercept)`) # remove Intercept column

# save
write.table(df_lasso, file.path(output_path, paste0(batch_num,".txt")), 
            quote=FALSE, sep="\t", row.names=TRUE)
