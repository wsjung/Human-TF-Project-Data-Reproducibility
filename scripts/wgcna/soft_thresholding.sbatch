#!/bin/bash
#SBATCH --job-name="wgcna_soft_threshold"
#SBATCH --output=/scratch/mblab/d.p.ruskin/human/logs/wgcna_%A_%a.out
#SBATCH --error=/scratch/mblab/d.p.ruskin/human/logs/wgcna_%A_%a.err
#SBATCH --mem=8G
#SBATCH --array=1-20

# Load necessary modules
eval $(spack load --sh singularityce)

# Define paths
tissue='breast_mammary_tissue' # CHANGE BASED ON TISSUE WE'RE ANALYZING
input_file="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/${tissue}/gene_counts_vst_wgcna_genie3.csv"
output_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/wgcna_2024_08_08/${tissue}/soft_thresholding_intermediate_output"
mkdir -p ${output_dir}
final_output_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/wgcna_2024_08_08/${tissue}/soft_thresholding_final_output"
mkdir -p ${final_output_dir}
final_output_file="${final_output_dir}/soft_thresholding_results.csv"

# Beta value based on the SLURM_ARRAY_TASK_ID
beta_value=${SLURM_ARRAY_TASK_ID}

# Output file for each beta value in the intermediate output path
output_file="${output_dir}/soft_thresholding_beta_${beta_value}.csv"

# Run the R script
singularity exec --no-home \
    -B /scratch/mblab/d.p.ruskin/ \
    /scratch/mblab/d.p.ruskin/singularity/containers/rstudio-4.2.2_1.2.sif \
    Rscript /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/SCRIPTS/DTO_INPUT/wgcna_2024_08_08/soft_thresholding.R $input_file $output_file $beta_value

# After all 20 iterations as separate jobs are finished, this will concatenate the results into a final CSV file
if [ "$SLURM_ARRAY_TASK_ID" -eq 20 ]; then
    sleep 120  # Wait for a minute to ensure all jobs are finished
    awk FNR!=1 ${output_dir}/soft_thresholding_beta_*.csv > $final_output_file
    sort -k1 -n -t ',' -o $final_output_file $final_output_file
    echo "All results concatenated into $final_output_file"
fi

