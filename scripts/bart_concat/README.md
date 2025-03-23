Steps for running BART on concatenated data (tissue-agnostic)

1. Run './create_batches.sh'
3. Run 'sbatch bart.sbatch'
3. Make sure that there are the correct number of batch results in the output directory (should be 200 batch .txt files)
4. Run './concat_batches.sh'
