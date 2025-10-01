#!/bin/bash

#SBATCH -J filter_sig_qval_peaks_1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-799%198
#SBATCH -o ./logs/filter_sig_qval_peaks_1_%A_%a.out

source activate human_tf_project

input_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1/
output_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01/
filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt )

echo $filename

python3 -u filter_sig_qval_peaks_1.py \
    --input_path $input_path \
    --filename $filename \
    --output_path $output_path
