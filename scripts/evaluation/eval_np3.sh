#!/bin/bash

#
#   This script runs all evaluation metrics for the NP3 network
#

tissue=$1
echo "EVALUATING TISSUE: ${tissue}"

p_src_code=/scratch/mblab/jungw/NET-evaluation/
p_src_res=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/res/
p_src_data=/scratch/mblab/jungw/human_TF_project/data/NP3_INPUT/tissue/${tissue}/model12/
p_out_dir=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/evaluation/
p_singularity_img=/scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif
p_singularity_bindpath=/scratch/mblab/jungw/

# BINDING
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
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# GO
#  requires that the input are stable IDs
input_network=${p_src_res}/net_np3.tsv
stable_network=${p_src_res}/net_np3_stableID.tsv
./stabilize_ensg_id.sh $input_network $stable_network
${p_src_code}go \
    --p_in_net "${stable_network}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go/eval_np3_go.tsv" \
    --p_out_dir "${p_out_dir}/go/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath


# go-directness command
#   requires that the binding labels are stableIDs to match stable ID results from GO evaluation
${p_src_code}go_directness \
    --p_in_dir_go "${p_out_dir}/go/" \
    --nbr_top_edges 50 \
    --nbr_edges_per_threshold 5 \
    --p_binding_event "${p_src_data}binding_labels_stableID.txt" \
    --flag_debug "ON" \
    --p_out_eval "${p_out_dir}/eval_np3_go_directness.tsv" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath


# ppi command
input_network=${p_src_res}/net_np3.tsv
stable_network=${p_src_res}/net_np3_stableID.tsv
./stabilize_ensg_id.sh $input_network $stable_network
${p_src_code}ppi \
    --p_in_net "${stable_network}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_np3_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath
