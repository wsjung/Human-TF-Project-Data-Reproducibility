#!/bin/bash

#SBATCH --mem=3G
#SBATCH --output=./logs/"find_remap_tfs_%A.txt"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu

conda_env_dir=$1
pipeline_dir=$2
scripts_dir=$3
remap_dir=$4

# Load conda environment
eval $(spack load --sh miniconda3)
source activate ${conda_env_dir}

python3 -u ${scripts_dir}find_remap_tfs.py \
	--master_path ${pipeline_dir}xgboost_master.tsv \
	--remap_path ${remap_dir} \
	--output_path ${pipeline_dir}

# Make 50 batches (if you want 100 batches, switch to "l/100 -a 3", then change number of array jobs in "eQTL_analysis_main_pyranges.sh")
mkdir -p ${pipeline_dir}remap_batches
split -n l/50 -a 2 --numeric-suffixes=1 --additional-suffix=.txt ${pipeline_dir}all_tfs.txt ${pipeline_dir}remap_batches/batch

echo 'batches saved to:'
echo ${pipeline_dir}remap_batches
