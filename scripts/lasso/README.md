# LASSO regression
==================

Scripts for using LASSO to predict gene expression using expression of all TFs.
Results from this regession are used to construct tissue-specific LASSO and
tissue-aggregate LASSO features for the XGBoost model.

## Input
--------

This requires you have an input directory with a subdirectory for each of the
tissues. Within each tissue subdirectory, there should be `np3_tpm.txt` file and
a `np3_tpm_tfs.txt` file.

* The `np3_tpm.txt` file is a tab-delimited file of gene TPM values with gene IDs
as rownames and sample IDs as column names.
* The `np3_tpm_tfs.txt` file is a subset of `np3_tpm.txt` with only the TFs.

Note that the input TPM data have been filtered for lowly expressed genes, as
per `../expression_filtering`.


## Steps
--------

1. Prepare batches of genes
    * `create_lasso_batches_1.py`
    * By preparing batches, we can run LASSO more quickly in parallel.
    * This script creates a set of .txt files listing which genes are grouped
      together as a batch.

2. Run LASSO
    * `lasso_2.R`
    * Runs the LASSO model for each gene in the specified batch.

3. Concatenate batch results to a single matrix
    * `lasso_concat_batch_results_3.R`

## Output
---------

For each tissue, the pipeline produces a `LASSO.txt` file, as well as one
`LASSO.txt` file for the tissue-aggregate version. Each file is tab-delimited,
with target genes as rows and TFs as columns. The entry in each cell corresponds
to the coefficient estimated by the LASSO regression model.

