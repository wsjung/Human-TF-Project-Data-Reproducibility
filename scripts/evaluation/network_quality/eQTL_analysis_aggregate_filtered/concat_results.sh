#!/bin/bash

#SBATCH --array=0-35
#SBATCH --mem=2G
#SBATCH --output=./logs/concat_results_%A_%a.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu

conda_env_dir=$1
pipeline_dir=$2
scripts_dir=$3
variant_gene_pairs_dir=$4
variant_type=$5

# Load conda environment
eval $(spack load --sh miniconda3)
source activate ${conda_env_dir}

python3 -u ${scripts_dir}concat_results.py \
	--master_path ${pipeline_dir}batch_outputs_per_eQTL_tissue \
	--variant_gene_pairs_path ${variant_gene_pairs_dir} \
	--variant_file_number ${SLURM_ARRAY_TASK_ID} \
	--variant_type ${variant_type} \
	--output_path ${pipeline_dir}
