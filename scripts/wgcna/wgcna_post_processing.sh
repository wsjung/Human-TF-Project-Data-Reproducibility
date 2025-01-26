#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH --array=0-47
#SBATCH -o /scratch/mblab/d.p.ruskin/human/logs/wgcna_post_processing_%A_%a.out

# Do a job for each tissue in the directory
TISSUE_DIRS=(/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/wgcna_no_design_2024_09_11/*)
TISSUE_DIR=${TISSUE_DIRS[$SLURM_ARRAY_TASK_ID]}
echo ${TISSUE_DIR}

# Load necessary modules
eval $(spack load --sh singularityce)

# Run the R script
singularity exec --no-home \
    -B /scratch/mblab/d.p.ruskin/ \
    /scratch/mblab/d.p.ruskin/singularity/containers/rstudio-4.2.2_1.2.sif \
    Rscript /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/SCRIPTS/DTO_INPUT/wgcna_2024_08_08/wgcna_post_processing.R \
    ${TISSUE_DIR}/wgcna_output \
    /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/REMAP_HGNC_ENSG_mapping.txt

