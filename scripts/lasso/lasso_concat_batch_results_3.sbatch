#!/bin/bash

#SBATCH -J lasso_concat_batch_results
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-37
#SBATCH -o ./logs/lasso_concat_batch_results_%A_%a.out

export R_LIBS_SITE=/scratch/mblab/jungw/renv/4.4/
export R_LIBS=/scratch/mblab/jungw/renv/4.4/
export TMP=/scratch/mblab/jungw/R_TMP/
export TMPDIR=/scratch/mblab/jungw/R_TMP/
eval $(spack load --sh r@4.4.0)

tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p ./tissues.txt )

# REMAP TF SET
results_path=/scratch/mblab/jungw/human_TF_project/data/GTEx_RNAseq/GTEx_processing_gene_counts_tpm_v29_output_lasso/${tissue}/

Rscript lasso_concat_batch_results_3.R \
    $results_path
