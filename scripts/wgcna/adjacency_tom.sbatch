#!/bin/bash
#SBATCH --job-name="whole_blood_adjacency_tom"
#SBATCH --output=/scratch/mblab/d.p.ruskin/human/logs/adjacency_tom_%A.out
#SBATCH --error=/scratch/mblab/d.p.ruskin/human/logs/adjacency_tom_%A.err
#SBATCH --mem=20G

# Load necessary modules
eval $(spack load --sh singularityce)

# All paths
tissue='whole_blood'

input_file="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/${tissue}/no_design_gene_counts_vst_wgcna_genie3.csv"
output_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/wgcna_no_design_2024_09_11/${tissue}/wgcna_output"
mkdir -p ${output_dir}
adjacency_output_file="${output_dir}/adjacency_matrix.csv"
tom_output_file="${output_dir}/tom.csv"
diss_tom_output_file="${output_dir}/diss_tom.csv"

best_beta=8 # Set this to the best beta value obtained from previous script

# Run the R script
singularity exec --no-home \
    -B /scratch/mblab/d.p.ruskin/ \
    /scratch/mblab/d.p.ruskin/singularity/containers/rstudio-4.2.2_1.2.sif \
    Rscript /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/SCRIPTS/DTO_INPUT/wgcna_2024_08_08/adjacency_tom.R $input_file $adjacency_output_file $tom_output_file $diss_tom_output_file $best_beta

