#!/bin/bash

#SBATCH -J create_cv_folds_4
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-36#1-36
#SBATCH -o ./logs/create_cv_folds_4_%A_%a.out

source activate human_tf_project

tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/tissues.txt )
base_path=/scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/
binding_threshold_10perc=10
model="metanet"

python3 -u create_cv_folds_4.py \
    --tissue $tissue \
    --base_path $base_path \
    --model $model \
    --num_CV_folds 10 \
    --binding_threshold $binding_threshold_10perc
