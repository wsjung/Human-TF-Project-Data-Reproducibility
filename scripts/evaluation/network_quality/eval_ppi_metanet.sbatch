#!/bin/bash

#SBATCH -J eval_ppi_metanet
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-36
#SBATCH -o ./logs/eval_ppi_metanet_%A_%a.out

source activate metanet

tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p ./tissues.txt )

PPI_db=./NET-evaluation/metadata/human/9606.STRING.gencodev29_ENSG_700.txt
top_TF_pairs=100
pairs_per_threshold=10
genes_per_TF=250

base_dir=/scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/results/
data_dir="${base_dir}/${tissue}/METANet/${model}/"
input_net="${data_dir}/${tissue}_${model}_stableID.tsv"
output_path="${data_dir}/evaluation/eval_ppi_${model}.tsv"

python3 -u eval_ppi.py \
    --input_net $input_net \
    --PPI_db $PPI_db \
    --top_TF_pairs $top_TF_pairs \
    --pairs_per_threshold $pairs_per_threshold \
    --genes_per_TF $genes_per_TF \
    --output_path $output_path

echo "DONE"
