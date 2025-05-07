#!/bin/bash

#SBATCH --array=1-50
#SBATCH --mem=3G
#SBATCH --output=./logs/"eQTL_analysis_main_pyranges_%A_%a.txt"
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu

variant_file_number=$1
conda_env_dir=$2
pipeline_dir=$3
scripts_dir=$4
aggregate_file=$5
remap_dir=$6
variant_gene_pairs_dir=$7
variant_type=$8
gencode_ensg_mappings=$9

# Load conda environment
eval $(spack load --sh miniconda3)
source activate ${conda_env_dir}

batch_dir=${pipeline_dir}remap_batches/
batch_file=$(ls $batch_dir | sed -n ${SLURM_ARRAY_TASK_ID}p)

python3 -u ${scripts_dir}eQTL_analysis_main_pyranges.py \
	--master_path ${aggregate_file} \
	--remap_batch ${batch_dir}${batch_file} \
	--remap_path ${remap_dir} \
	--variant_gene_pairs_path ${variant_gene_pairs_dir} \
	--variant_file_number ${variant_file_number} \
	--variant_type ${variant_type} \
	--gencode_ensg_mapping_path ${gencode_ensg_mappings} \
	--output_path ${pipeline_dir}
