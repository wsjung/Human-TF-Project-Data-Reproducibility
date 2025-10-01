#!/bin/bash

#SBATCH -J lasso
#SBATCH --mem-per-cpu=12G
#SBATCH --array=3602-3800
#SBATCH -o ./logs/lasso_%A_%a.out

# TISSUE-SPECIFIC #
tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/scripts/xgboost_all_gtex_tissues/lasso_gtex/tissuesXbatches.txt | cut -f1 )
batch_num=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/scripts/xgboost_all_gtex_tissues/lasso_gtex/tissuesXbatches.txt | cut -f2 )

# REMAP TF SET
input_data_path=/scratch/mblab/jungw/human_TF_project/data/GTEx_RNAseq/GTEx_processing_gene_counts_tpm_v29_output/${tissue}/
output_data_path=/scratch/mblab/jungw/human_TF_project/data/GTEx_RNAseq/GTEx_processing_gene_counts_tpm_v29_output_lasso/${tissue}/batches/

export R_LIBS_SITE=/scratch/mblab/jungw/renv/4.4/
export R_LIBS=/scratch/mblab/jungw/renv/4.4/
export TMP=/scratch/mblab/jungw/R_TMP/
export TMPDIR=/scratch/mblab/jungw/R_TMP/
eval $(spack load --sh r@4.4.0)

Rscript lasso_2.R \
    $input_data_path \
    $batch_num \
    $output_data_path
