#!/bin/bash

#SBATCH -J max_remap_peak_sum_scores
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-799#1-799%198
#SBATCH -o ./logs/max_remap_peak_sum_scores_%A_%a.out

source activate human_tf_project


remap_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt )
input_path="/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping/${remap_file}"
output_path="/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_max_sum/${remap_file}"

python3 -u max_remap_peak_sum_scores_3.py \
    --input_path $input_path \
    --output_path $output_path

