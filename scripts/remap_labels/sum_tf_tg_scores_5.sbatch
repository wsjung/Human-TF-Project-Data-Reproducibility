#!/bin/bash

#SBATCH -J sum_tf_tg_scores
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-799%198
#SBATCH -o ./logs/sum_tf_tg_scores_%A_%a.out

source activate human_data


remap_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt )
input_file="/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_max_sum_annotated/${remap_file}"
output_file="/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_max_sum_annotated_TF_TG_scores/${remap_file}"
gencode=/scratch/mblab/jungw/human_data/gencode_v29/gencode.v29.protein_coding_genes_aliases.tsv

python3 -u sum_tf_tg_scores.py \
    --input_file $input_file \
    --output_file $output_file \
    --gencode $gencode

