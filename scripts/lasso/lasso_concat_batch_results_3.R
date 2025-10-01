library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
results_path <- args[1]
batches_path <- file.path(results_path, "batches")

# 1. list of batch result file paths
list_batch_filepaths <- list.files(
    batches_path,
    pattern="^batch\\d{3}\\.txt$",
    full.names=TRUE
)
print(paste0("Num batch files: ", length(list_batch_filepaths)))

# 2. collect batch results in a list
print("collecting batch results")
list_batch_dfs <- lapply(list_batch_filepaths, function(filepath) {
    df <- read.table(filepath, sep="\t")
    return(df)
})

# 3. concatenate dfs
print("concatenating")
df_final <- dplyr::bind_rows(list_batch_dfs)

# 4. write
output_path <- file.path(results_path,"LASSO.txt")
print(paste0("writing to ", output_path))
df_final %>% write.table(output_path, sep="\t", row.names=TRUE, quote=FALSE)

print("DONE")
