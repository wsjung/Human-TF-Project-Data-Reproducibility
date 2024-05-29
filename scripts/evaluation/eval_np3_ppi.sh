#!/bin/bash

#
#   This script runs the NP3 PPI evaluation for the NP3 network and all input features
#

tissue=$1
echo "EVALUATING TISSUE: ${tissue}"

p_src_code=/scratch/mblab/jungw/NET-evaluation/
p_src_res=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/res/
p_src_data=/scratch/mblab/jungw/human_TF_project/data/NP3_INPUT/tissue/${tissue}/model12/
p_out_dir=/scratch/mblab/jungw/human_TF_project/data/NP3_results/${tissue}/evaluation/
p_singularity_img=/scratch/mblab/jungw/NET-evaluation/singularity/s_neteval.sif
p_singularity_bindpath=/scratch/mblab/jungw/

# NP3 network
# convert NP3 network to stable gene IDs
input_network_np3=${p_src_res}/net_np3.tsv
stable_network_np3=${p_src_res}/net_np3_stableID.tsv
if [ ! -f ${stable_network_np3} ]; then
    ./stabilize_ensg_id.sh $input_network_np3 $stable_network_np3
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_np3}" \
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

# LASSO network
# convert NP3 network to stable gene IDs
input_network_lasso=${p_src_res}/tmp_combine/network_construction/net_lasso.tsv
stable_network_lasso=${p_src_res}/tmp_combine/network_construction/net_lasso_stableID.tsv
if [ ! -f ${stable_network_lasso} ]; then
    ./stabilize_ensg_id.sh $input_network_lasso $stable_network_lasso
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_lasso}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_lasso_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# BART network
# convert NP3 network to stable gene IDs
input_network_bart=${p_src_res}/tmp_combine/network_construction/net_bart.tsv
stable_network_bart=${p_src_res}/tmp_combine/network_construction/net_bart_stableID.tsv
if [ ! -f ${stable_network_bart} ]; then
    ./stabilize_ensg_id.sh $input_network_bart $stable_network_bart
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_bart}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_bart_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# DE network
# convert NP3 network to stable gene IDs
input_network_de=${p_src_res}/tmp_combine/network_construction/net_de.tsv
stable_network_de=${p_src_res}/tmp_combine/network_construction/net_de_stableID.tsv
if [ ! -f ${stable_network_de} ]; then
    ./stabilize_ensg_id.sh $input_network_de $stable_network_de
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_de}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_de_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# MOTIF(EA) network
# convert NP3 network to stable gene IDs
input_network_motifea=${p_src_res}/tmp_combine/network_construction/net_motifea.tsv
stable_network_motifea=${p_src_res}/tmp_combine/network_construction/net_motifea_stableID.tsv
if [ ! -f ${stable_network_motifea} ]; then
    ./stabilize_ensg_id.sh $input_network_motifea $stable_network_motifea
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_motifea}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_motifea_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# MOTIF(f5) network
# convert NP3 network to stable gene IDs
input_network_motiff5=${p_src_res}/tmp_combine/network_construction/net_motiff5.tsv
stable_network_motiff5=${p_src_res}/tmp_combine/network_construction/net_motiff5_stableID.tsv
if [ ! -f ${stable_network_motiff5} ]; then
    ./stabilize_ensg_id.sh $input_network_motiff5 $stable_network_motiff5
fi
${p_src_code}ppi \
    --p_in_net "${stable_network_motiff5}" \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/eval_motiff5_ppi.tsv" \
    --nbr_top_edges 100 \
    --nbr_edges_per_threshold 10 \
    --threshold 25 \
    --p_STRING_db "${p_src_code}metadata/human/STRINGdb_PPI_ensg_high_confidence_700.txt" \
    --STRING_confidence 700 \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath
