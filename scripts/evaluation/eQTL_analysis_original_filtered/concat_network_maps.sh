#!/bin/bash

#SBATCH --mem=8G
#SBATCH --output=./logs/concat_network_maps_%A.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu

conda_env_dir=$1
pipeline_dir=$2
scripts_dir=$3
human_network_files_dir=$4

# Load conda environment
eval $(spack load --sh miniconda3)
source activate ${conda_env_dir}

python3 -u ${scripts_dir}concat_network_maps.py \
	--human_network_files_path ${human_network_files_dir} \
	--output_path ${pipeline_dir}

