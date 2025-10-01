#!/bin/bash

#
#   This script runs the NP3 GO-directness evaluation for the NP3 network and all its input feature networks
#

tissue=$1
echo "EVALUATING TISSUE: ${tissue}"

p_src_code=/scratch/mblab/jungw/NET-evaluation/
p_src_res=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/res/
p_src_data=/scratch/mblab/jungw/human_TF_project/data/NP3_INPUT/tissue/${tissue}/model12/
p_out_dir=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/evaluation/
p_singularity_bindpath=/scratch/mblab/jungw/
p_singularity_img=/scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif

# NP3 network
# convert NP3 network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_np3/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_debug "ON" \
    --p_out_eval "${p_out_dir}/go_np3/eval_np3_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# LASSO
# convert LASSO network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_lasso/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_debug "ON" \
    --p_out_eval "${p_out_dir}/go_lasso/eval_lasso_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# BART 
# convert BART network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_bart/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_debug "ON" \
    --p_out_eval "${p_out_dir}/go_bart/eval_bart_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# DE
# convert DE network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_de/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_debug "ON" \
    --p_out_eval "${p_out_dir}/go_de/eval_de_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath


# MOTIF (EA)
# convert MOTIF (EA) network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_motifea/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_motifeabug "ON" \
    --p_out_eval "${p_out_dir}/go_motifea/eval_motifea_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath


# MOTIF (F5)
# convert MOTIF (F5) network to stable gene IDs
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go_motiff5/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_motiff5bug "ON" \
    --p_out_eval "${p_out_dir}/go_motiff5/eval_motiff5_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath
