#!/bin/bash

#SBATCH -J eval_binding_metanet
#SBATCH --mem=3G
#SBATCH --array=1-36
#SBATCH -o ./logs/eval_binding_metanet_%A_%a.out

source activate metanet

tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p ./tissues.txt )
echo $tissue

# paths
base_path="./data/"
model="METANet"
submodel="xgboost"

network_path="${base_path}/results/${tissue}/${model}/${submodel}/${tissue}_${submodel}.tsv"
label_path="${base_path}/input_data/${tissue}/${model}/features/binding_labels.txt"
output_path="${base_path}/results/${tissue}/${model}/${submodel}/evaluation/eval_binding_${submodel}.tsv"

max_avg_targets=250

python3 -u eval_binding.py \
    --network $network_path \
    --label $label_path \
    --max_avg_targets $max_avg_targets \
    --output $output_path

echo "DONE"
