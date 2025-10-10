#!/bin/bash

#SBATCH -J eval_go250_metanet
#SBATCH --array=1-36
#SBATCH --mem-per-cpu=2G
#SBATCH -o ./logs/eval_go250_metanet_%A_%a.out

### update this ###
eval_binding_script=./eval_go250_and_directness_fixed_penalty_corrected_pval.sh

### update this ###
tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p ./tissues.txt )

### update this ###
net_eval_dir=$HOME/NET-evaluation/
singularity_bindpath=$HOME
singularity_img="${net_eval_dir}/singularity/s_neteval.sif"

basepath=$HOME/human_TF_project/data/xgboost_all_gtex_tissues/

results_dir="${basepath}/results/${tissue}/METANEt/"
features_dir="${basepath}/input_data/${tissue}/METANet/features/"

# 10% threshold #
${eval_binding_script} $tissue $results_dir $features_dir $net_eval_dir $singularity_bindpath $singularity_img $stabilize_gene
