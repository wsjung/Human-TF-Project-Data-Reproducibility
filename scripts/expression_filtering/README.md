# Filtering of GTEx expression data

* The scripts in this directory serve to remove lowly expressed genes.
* The filtering is carried out using gene read counts and the filtered results
  are used to select genes in the TPM dataset.
    * Genes with <= 3 CPM in >= 98.5% of samples were removed.
* `high_expr_filter_tpm.R` filters tissue-specific expression data.
* `concat_high_expr_filter_tpm.R` filters tissue-aggregate expression data
  (samples across all tissues).
