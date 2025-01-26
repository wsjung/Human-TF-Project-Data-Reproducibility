#!/bin/bash

#SBATCH --array=1-199
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/mblab/d.p.ruskin/human/logs/genie3_%A_%a.txt
#SBATCH -e /scratch/mblab/d.p.ruskin/human/logs/genie3_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu
#SBATCH --job-name=genie3_vagina

# Load conda environment
eval $(spack load --sh miniconda3)
source activate /scratch/mblab/d.p.ruskin/conda/conda_environments/human_network

tissue="vagina"

# Path to batch files
batch_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/${tissue}/genie3_batches"
output_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/${tissue}/genie3_output"
mkdir -p ${output_dir}

# Get batch file corresponding to the current task
batch_file=$(ls $batch_dir | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Run the Python script for the current batch file
python3 -u /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/SCRIPTS/DTO_INPUT/genie3_2024_08_08/genie3.py --input_file /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/${tissue}/no_design_gene_counts_vst_wgcna_genie3.csv --batch_file $batch_dir/$batch_file --output_dir ${output_dir}
