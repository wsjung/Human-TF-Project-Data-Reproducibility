#!/bin/bash

#SBATCH -J create_lasso_batches.py
#SBATCH --mem=8G
#SBATCH --array=1-36#1-1
#SBATCH -o ./logs/create_lasso_batches_%A_%a.out

source activate human_tf_project

# TF set based on LAMBERT #
input_path=/scratch/mblab/jungw/human_TF_project/data/GTEx_RNAseq/GTEx_processing_gene_counts_tpm_v29_output_all_lambert_TFs/

# tissue specific #
tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/tissues.txt )

## tissue-aggregate #
#tissue="tpm_concat_output"
#    # -> use --aggregate flag

python3 -u create_lasso_batches_1.py \
    --input_path $input_path \
    --tissue $tissue
