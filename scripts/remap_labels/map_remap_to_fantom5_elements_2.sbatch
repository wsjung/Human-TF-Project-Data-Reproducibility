#!/bin/bash

#SBATCH -J map_remap_to_fantom5_elements
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-799%198
#SBATCH -o ./logs/map_remap_to_fantom5_elements_%A_%a.out

source activate human_tf_project

fantom5_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/fantom5_enhancer_promoter_overlapped_regions_only.bed
remap_input_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01/
remap_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt )
remap_path="${remap_input_path}/${remap_file}"
output_path="/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping/${remap_file}"


python3 map_remap_to_fantom5_elements_2.py \
    --fantom5_path $fantom5_path \
    --remap_path $remap_path \
    --output_path $output_path
