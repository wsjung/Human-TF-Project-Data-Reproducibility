library(dplyr)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
tissue <- args[1]
base_path <- args[2]
binding_threshold <- as.numeric(args[3])


plot_dir <- file.path(base_path, "plots", tissue, "feature_transformation_metanet")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
all_cols <- c("MARBACH","BART_GTEX","BART_ALL_TISSUES","LASSO_GTEX","LASSO_ALL_TISSUES")
print(plot_dir)

# load data
print("loading data")
df_path <- file.path(base_path, "input_data", tissue, paste0(binding_threshold, ".0_threshold_metanet"), paste0(tissue, "_", binding_threshold, ".0_threshold_metanet.txt"))
df = read.table(df_path, sep="\t", header=TRUE)
print(head(df))


### plot original variables
print("plot original")
df_long <- df %>%
  pivot_longer(cols = all_of(all_cols),
               names_to = "original",
               values_to = "value")
ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~original, scales = "free", nrow=4, ncol=3) +
  labs(title = "Original feature distributions",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "original_features.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~original, scales = "free", nrow=4, ncol=3) +
  labs(title = "Original feature distributions (log10 y-axis)",
       x = "Value",
       y = "Frequency (log10)")
ggsave(file.path(plot_dir, "original_features_log10y.pdf"), device="pdf")


### MARBACH
print("marbach")
df$MARBACH_rank <- rank(desc(df$MARBACH), ties.method = 'max')
df$MARBACH_nlogrank <- log10(df$MARBACH_rank) * -1
df$MARBACH_nlogrank_minmax <- (df$MARBACH_nlogrank - min(df$MARBACH_nlogrank)) / (max(df$MARBACH_nlogrank) - min(df$MARBACH_nlogrank))

df_long <- df %>%
  pivot_longer(cols = starts_with("MARBACH"),
               names_to = "transformation",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "MARBACH feature transformation",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "MARBACH_nlogrank_transformations.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "MARBACH feature transformation (log10 y)",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "MARBACH_nlogrank_transformations_log10y.pdf"), device="pdf")


### BART_GTEX
print("BART TPM")
df$BART_GTEX_rank <- rank(desc(df$BART_GTEX), ties.method = 'max')
df$BART_GTEX_nlogrank <- log10(df$BART_GTEX_rank) * -1
df$BART_GTEX_nlogrank_minmax <- (df$BART_GTEX_nlogrank - min(df$BART_GTEX_nlogrank)) / (max(df$BART_GTEX_nlogrank) - min(df$BART_GTEX_nlogrank))

df_long <- df %>%
  pivot_longer(cols = starts_with("BART_GTEX"),
               names_to = "transformation",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "BART_GTEX feature transformation",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "BART_GTEX_nlogrank_transformations.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "BART_GTEX feature transformation (log10 y)",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "BART_GTEX_nlogrank_transformations_log10y.pdf"), device="pdf")


### LASSO_GTEX
print("LASSO TPM")
df$LASSO_GTEX_rank <- rank(desc(df$LASSO_GTEX), ties.method = 'max')
df$LASSO_GTEX_nlogrank <- log10(df$LASSO_GTEX_rank) * -1
df$LASSO_GTEX_nlogrank_minmax <- (df$LASSO_GTEX_nlogrank - min(df$LASSO_GTEX_nlogrank)) / (max(df$LASSO_GTEX_nlogrank) - min(df$LASSO_GTEX_nlogrank))

df_long <- df %>%
  pivot_longer(cols = starts_with("LASSO_GTEX"),
               names_to = "transformation",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "LASSO_GTEX feature transformation",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "LASSO_GTEX_nlogrank_transformations.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "LASSO_GTEX feature transformation (log10 y)",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "LASSO_GTEX_nlogrank_transformations_log10y.pdf"), device="pdf")



### BART_ALL_TISSUES
print("BART ALL GTEX")
df$BART_ALL_TISSUES_rank <- rank(desc(df$BART_ALL_TISSUES), ties.method = 'max')
df$BART_ALL_TISSUES_nlogrank <- log10(df$BART_ALL_TISSUES_rank) * -1
df$BART_ALL_TISSUES_nlogrank_minmax <- (df$BART_ALL_TISSUES_nlogrank - min(df$BART_ALL_TISSUES_nlogrank)) / (max(df$BART_ALL_TISSUES_nlogrank) - min(df$BART_ALL_TISSUES_nlogrank))

df_long <- df %>%
  pivot_longer(cols = starts_with("BART_ALL_TISSUES"),
               names_to = "transformation",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "BART_ALL_TISSUES feature transformation",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "BART_ALL_TISSUES_nlogrank_transformations.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "BART_ALL_TISSUES feature transformation (log10 y)",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "BART_ALL_TISSUES_nlogrank_transformations_log10y.pdf"), device="pdf")


### LASSO_ALL_TISSUES
print("LASSO ALL GTEX")
df$LASSO_ALL_TISSUES_rank <- rank(desc(df$LASSO_ALL_TISSUES), ties.method = 'max')
df$LASSO_ALL_TISSUES_nlogrank <- log10(df$LASSO_ALL_TISSUES_rank) * -1
df$LASSO_ALL_TISSUES_nlogrank_minmax <- (df$LASSO_ALL_TISSUES_nlogrank - min(df$LASSO_ALL_TISSUES_nlogrank)) / (max(df$LASSO_ALL_TISSUES_nlogrank) - min(df$LASSO_ALL_TISSUES_nlogrank))

df_long <- df %>%
  pivot_longer(cols = starts_with("LASSO_ALL_TISSUES"),
               names_to = "transformation",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "LASSO_ALL_TISSUES feature transformation",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "LASSO_ALL_TISSUES_nlogrank_transformations.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~ transformation, scales = "free_x") +
  labs(title = "LASSO_ALL_TISSUES feature transformation (log10 y)",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "LASSO_ALL_TISSUES_nlogrank_transformations_log10y.pdf"), device="pdf")



### nlogrank transformed features
print("FINAL")
df_nlogrank_minmax <- df[,c("TF","GENE","LABEL")]
df_nlogrank_minmax$MARBACH <- df$MARBACH_nlogrank_minmax
df_nlogrank_minmax$BART_GTEX <- df$BART_GTEX_nlogrank_minmax
df_nlogrank_minmax$BART_ALL_TISSUES <- df$BART_ALL_TISSUES_nlogrank_minmax
df_nlogrank_minmax$LASSO_GTEX <- df$LASSO_GTEX_nlogrank_minmax
df_nlogrank_minmax$LASSO_ALL_TISSUES <- df$LASSO_ALL_TISSUES_nlogrank_minmax


### plot final features
df_long <- df_nlogrank_minmax %>%
  pivot_longer(cols = all_cols,
               names_to = "features",
               values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) +
  facet_wrap(~features, scales = "free", nrow=4, ncol=3) +
  labs(title = "Transformed feature distributions",
       x = "Value",
       y = "Frequency")
ggsave(file.path(plot_dir, "transformed_features_nlogrank_minmax.pdf"), device="pdf")

ggplot(df_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black", bins = 30) + scale_y_continuous(trans='log10') +
  facet_wrap(~features, scales = "free", nrow=4, ncol=3) +
  labs(title = "Transformed feature distributions (log10 y-axis)",
       x = "Value",
       y = "Frequency (log10)")
ggsave(file.path(plot_dir, "transformed_features_nlogrank_minmax_log10y.pdf"), device="pdf")

## save file
output_filepath = file.path(base_path, "input_data", tissue, paste0(binding_threshold, ".0_threshold_metanet"), paste0(tissue, "_", binding_threshold, ".0_threshold_metanet_nlogrank_minmax.txt"))
write.table(df_nlogrank_minmax, file=output_filepath, sep="\t", row.names=FALSE, quote=FALSE)


print("DONE")
