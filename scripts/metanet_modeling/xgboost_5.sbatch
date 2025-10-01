#!/bin/bash

#SBATCH -J xgboost
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=5
#SBATCH --array=1-36
#SBATCH -o ./logs/xgboost_5_%A_%a.out
#SBATCH --mail-type=END,FAIL

source activate human_tf_project

###########################################
### update this to your list of tissues ###
###########################################
tissues_file=/scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/tissues.txt

tissue=$( sed -n ${SLURM_ARRAY_TASK_ID}p $tissues_file )

########################
### UPDATE BASE PATH ###
########################
base_path=/scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/

n_cv_folds=10
seed=42

binding_threshold_10perc=10

###################
### update this ###
###################
binding_threshold=$binding_threshold_10perc

model="metanet"
output_subdir="xgboost"    # testing

python3 -u xgboost_5.py \
    --tissue $tissue \
    --base_path $base_path \
    --model $model \
    --n_cv_folds $n_cv_folds \
    --seed $seed \
    --binding_threshold $binding_threshold \
    --output_subdir $output_subdir
