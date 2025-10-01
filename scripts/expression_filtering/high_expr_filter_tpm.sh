#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-53
#SBATCH -o /scratch/mblab/d.p.ruskin/human/logs/high_expr_filter_tpm_%A_%a.out

# Do a job for each tissue in the directory
TISSUE_TPM_DIRS=(/scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/GTEx_processing_gene_counts_tpm_v29_output/*)
TISSUE_TPM_DIR=${TISSUE_TPM_DIRS[$SLURM_ARRAY_TASK_ID]}
echo ${TISSUE_TPM_DIR}

TISSUE_GENE_COUNTS_DIRS=(/scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/*)   
TISSUE_GENE_COUNTS_DIR=${TISSUE_GENE_COUNTS_DIRS[$SLURM_ARRAY_TASK_ID]}
echo ${TISSUE_GENE_COUNTS_DIR}

# Load necessary modules
eval $(spack load --sh singularityce)

# Run the R script
singularity exec --no-home \
    -B /scratch/mblab/d.p.ruskin/ \
    /scratch/mblab/d.p.ruskin/singularity/containers/rstudio-4.2.2_1.2.sif \
    Rscript /scratch/mblab/d.p.ruskin/human/scripts/protein_coding_high_expr_filter/high_expr_filter_tpm.R \
    ${TISSUE_TPM_DIR} \
    ${TISSUE_GENE_COUNTS_DIR} \
    /scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/Lambert_TFs.csv

