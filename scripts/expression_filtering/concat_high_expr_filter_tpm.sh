#!/bin/bash

#SBATCH --mem-per-cpu=32G
#SBATCH -o /scratch/mblab/d.p.ruskin/human/logs/concat_high_expr_filter_tpm_%A.out

# Load necessary modules
eval $(spack load --sh singularityce)

# Run the R script
singularity exec --no-home \
    -B /scratch/mblab/d.p.ruskin/ \
    /scratch/mblab/d.p.ruskin/singularity/containers/rstudio-4.2.2_1.2.sif \
    Rscript /scratch/mblab/d.p.ruskin/human/scripts/protein_coding_high_expr_filter/concat_high_expr_filter_tpm.R \
    /scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/GTEx_processing_gene_counts_tpm_v29_output/ \
    /scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/ \
    /scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/Lambert_TFs.csv \
    /scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/tpm_concat_output

