#!/bin/bash

#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/mblab/d.p.ruskin/human/logs/genie3_post_processing_%A.txt
#SBATCH -e /scratch/mblab/d.p.ruskin/human/logs/genie3_post_processing_%A.err
#SBATCH --job-name=genie3_pp_whole_blood

# Load conda environment
eval $(spack load --sh miniconda3)
source activate /scratch/mblab/d.p.ruskin/conda/conda_environments/human_network

tissue="whole_blood"

input_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/${tissue}/"
output_dir="/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07_post_processing/${tissue}/genie3_post_processing/"
mkdir -p ${output_dir}

# Run a Python script which combines all the batch outputs
python3 -u /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/SCRIPTS/DTO_INPUT/genie3_2024_08_08/genie3_post_processing.py \
	--genie3_output ${input_dir}genie3_output \
	--genes_list ${input_dir}all_genes_newline_separated.txt \
	--tfs_list /scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/REMAP_HGNC_ENSG_mapping.txt \
	--output_dir ${output_dir} \
