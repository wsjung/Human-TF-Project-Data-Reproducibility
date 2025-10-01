#!/bin/bash

#
#   This script runs the NP3 GO evaluation for the NP3 network and all its input feature networks
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
input_network_np3=${p_src_res}/net_np3.tsv
stable_network_np3=${p_src_res}/net_np3_stableID.tsv
./stabilize_ensg_id.sh $input_network_np3 $stable_network_np3
${p_src_code}go \
    --p_in_net "${stable_network_np3}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go_np3/eval_np3_go.tsv" \
    --p_out_dir "${p_out_dir}/go_np3/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# LASSO
# convert LASSO network to stable gene IDs
input_network_lasso=${p_src_res}/tmp_combine/network_construction/net_lasso.tsv
stable_network_lasso=${p_src_res}/tmp_combine/network_construction/net_lasso_stableID.tsv
./stabilize_ensg_id.sh $input_network_lasso $stable_network_lasso
${p_src_code}go \
    --p_in_net "${stable_network_lasso}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go_lasso/eval_lasso_go.tsv" \
    --p_out_dir "${p_out_dir}/go_lasso/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# BART 
# convert BART network to stable gene IDs
input_network_bart=${p_src_res}/tmp_combine/network_construction/net_bart.tsv
stable_network_bart=${p_src_res}/tmp_combine/network_construction/net_bart_stableID.tsv
./stabilize_ensg_id.sh $input_network_bart $stable_network_bart
${p_src_code}go \
    --p_in_net "${stable_network_bart}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go_bart/eval_bart_go.tsv" \
    --p_out_dir "${p_out_dir}/go_bart/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# DE
# convert DE network to stable gene IDs
input_network_de=${p_src_res}/tmp_combine/network_construction/net_de.tsv
stable_network_de=${p_src_res}/tmp_combine/network_construction/net_de_stableID.tsv
./stabilize_ensg_id.sh $input_network_de $stable_network_de
${p_src_code}go \
    --p_in_net "${stable_network_de}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go_de/eval_de_go.tsv" \
    --p_out_dir "${p_out_dir}/go_de/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# MOTIF (EA)
# convert MOTIF (EA) network to stable gene IDs
input_network_motifea=${p_src_res}/tmp_combine/network_construction/net_motifea.tsv
stable_network_motifea=${p_src_res}/tmp_combine/network_construction/net_motifea_stableID.tsv
./stabilize_ensg_id.sh $input_network_motifea $stable_network_motifea
${p_src_code}go \
    --p_in_net "${stable_network_motifea}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go_motifea/eval_motifea_go.tsv" \
    --p_out_dir "${p_out_dir}/go_motifea/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath

# MOTIF (F5)
# convert MOTIF (F5) network to stable gene IDs
input_network_motiff5=${p_src_res}/tmp_combine/network_construction/net_motiff5.tsv
stable_network_motiff5=${p_src_res}/tmp_combine/network_construction/net_motiff5_stableID.tsv
./stabilize_ensg_id.sh $input_network_motiff5 $stable_network_motiff5
${p_src_code}go \
    --p_in_net "${stable_network_motiff5}" \
    --p_gene_association "${p_src_code}metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${p_src_code}metadata/human/go-basic.obo" \
    --nbr_genes 15705 \
    --flag_slurm "ON" \
    --p_out_logs "${p_out_dir}logs/" \
    --data "${tissue}" \
    --p_out_eval "${p_out_dir}/go/eval_motiff5_go.tsv" \
    --p_out_dir "${p_out_dir}/go/" \
    --flag_singularity "ON" \
    --p_singularity_img $p_singularity_img \
    --p_singularity_bindpath $p_singularity_bindpath
