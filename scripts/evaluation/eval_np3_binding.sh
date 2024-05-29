#!/bin/bash

#
#   This script runs the binding support evaluation metric for the NP3 network and all input feature networks
#

tissue=$1
echo "EVALUATING TISSUE: ${tissue}"

p_src_code=/scratch/mblab/jungw/NET-evaluation/
p_src_res=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/res/
p_src_data=/scratch/mblab/jungw/human_TF_project/data/NP3_INPUT/tissue/${tissue}/model12/
p_out_dir=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/evaluation/
p_singularity_img=/scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif
p_singularity_bindpath=/scratch/mblab/jungw/

# FULL MODEL
${p_src_code}binding \
    --p_in_net "${p_src_res}/net_np3.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_np3_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img  $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# LASSO ONLY -- need to convert matrix from wide to long (already done in tmp_combine/network_construction)
${p_src_code}binding \
    --p_in_net "${p_src_res}/tmp_combine/network_construction/net_lasso.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_lasso_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img  $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# BART ONLY -- need to convert matrix from wide to long
${p_src_code}binding \
    --p_in_net "${p_src_res}/tmp_combine/network_construction//net_bart.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_bart_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img /scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif \
    --p_singularity_bindpath /scratch/mblab/jungw/

# DE ONLY -- need to convert matrix from wide to long
${p_src_code}binding \
    --p_in_net "${p_src_res}/tmp_combine/network_construction/net_de.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_de_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img /scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif \
    --p_singularity_bindpath /scratch/mblab/jungw/

# MOTIF (EA) ONLY -- need to convert matrix from wide to long
${p_src_code}binding \
    --p_in_net "${p_src_res}/tmp_combine/network_construction/net_motifea.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_motifea_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img  $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# MOTIF (F5) ONLY -- need to convert matrix from wide to long
${p_src_code}binding \
    --p_in_net "${p_src_res}/tmp_combine/network_construction/net_motiff5.tsv" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels.txt" \
    --p_out_eval "${p_out_dir}/eval_motiff5_binding.tsv" \
    --data "${tissue}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}/logs/" \
    --flag_singularity "ON" \
    --p_singularity_img  $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath
