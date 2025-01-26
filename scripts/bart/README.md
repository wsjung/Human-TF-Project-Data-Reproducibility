Steps for running BART on each GTEx tissue with filtered tpm data (example is liver)

1. Run './create_batches.sh liver'
2. Set tissue name to liver through 'tissue="liver"' line in bart.sbatch
3. Run 'sbatch bart.sbatch'
3. Make sure that there are the correct number of batch results in the output directory (should be 100 batch .txt files)
4. Run './concat_batches.sh liver'
