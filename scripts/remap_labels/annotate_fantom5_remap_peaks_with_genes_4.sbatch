#!/bin/bash

#SBATCH -J annotate_fantom5_remap_peaks_with_genes
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-799%198
#SBATCH -o ./logs/annotate_fantom5_remap_peaks_with_genes_%A_%a.out

source activate human_data


remap_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt )
tf_ensg=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_tfs.txt | cut -d . -f1-2)
input_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_max_sum/
output_path=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_max_sum_annotated/
gene_annotation=/scratch/mblab/jungw/human_TF_project/data/new_REMAP_binding_labels/fantom5_enhancer_promoter_overlapped_gene_annotation.bed

python3 -u annotate_fantom5_remap_peaks_with_genes.py \
    --input_path $input_path \
    --filename $remap_file \
    --TF_ensg $tf_ensg \
    --gene_annotation $gene_annotation \
    --output_path $output_path

