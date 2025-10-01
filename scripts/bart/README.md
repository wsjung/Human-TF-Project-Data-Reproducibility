# BART regression
==================

Scripts for using BART to infer TF-gene scores from expression data.
For each target gene, a BART model is trained to predict its expression level from the
experssion levels of all TFs.  Edge weights are then computed by virtually perturbing
each TF from its observed minimum to maximum (holding all other TFs at their
median) and taking the change in the model-predicted target expression. These
BART scores are later used as features in XGBoost.


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

Steps for running BART on each GTEx tissue with filtered tpm data (example is liver)

1. Run `./create_batches.sh liver`
2. Set tissue name to liver through 'tissue="liver"' line in bart.sbatch
3. Run `sbatch bart.sbatch`
3. Make sure that there are the correct number of batch results in the output directory (should be 100 batch .txt files)
4. Run `./concat_batches.sh liver`

## Output
---------

For each tissue, the pipeline produces a `${tissue}/bart_output/net_bart.tsv`
file.
    * A tab-delimited matrix with rows = TFs and columns = all target genes in
      the tissue. Each cell is the BART score defined above.
