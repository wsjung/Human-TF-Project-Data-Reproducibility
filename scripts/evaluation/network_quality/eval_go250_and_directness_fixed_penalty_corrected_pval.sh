#!/bin/bash
#
#   this script runs the GO and GO-directness evaluations
#

tissue=$1
results_dir=$2
features_dir=$3
net_eval_dir=$4
singularity_bindpath=$5
singularity_img=$6
stabilize_gene=$7

binding_label=$features_dir/binding_labels.txt
binding_label_stableid=$features_dir/binding_labels_stableID.txt

echo $tissue
echo $results_dir
echo $features_dir
echo $net_eval_dir
echo $singularity_bindpath
echo $singularity_img
echo $stabilize_gene
echo $binding_label
echo $binding_label_stableid
echo $evaluation_dir
echo $log_dir

###########
# XGBOOST #
###########
evaluation_dir="${results_dir}/xgboost/evaluation/"
log_dir="${results_dir}/xgboost/evaluation/logs/"
xgboost_fullid="${results_dir}/xgboost/${tissue}_xgboost.tsv"
xgboost_stableid="${results_dir}/xgboost/${tissue}_xgboost_stableID.tsv"

echo $xgboost_fullid
echo $xgboost_stableid

# GO
xgboost_jobid=$( ${net_eval_dir}/go_250 \
    --p_in_net ${xgboost_stableid} \
    --p_gene_association "${net_eval_dir}/metadata/human/gene_association.goa_human_ensembl.edit" \
    --p_gene_ontology "${net_eval_dir}/metadata/human/gene_ontology_edit.obo" \
    --p_out_dir "${evaluation_dir}/go250_xgboost_fixed_penalty_corrected_pval/" \
    --p_out_eval "${evaluation_dir}/go250_xgboost_fixed_penalty_corrected_pval/eval_go250_xgboost.tsv" \
    --p_out_logs "${log_dir}" \
    --nbr_genes 19927 \
    --data "go250_${tissue}_xgboost" \
    --flag_slurm "ON" \
    --flag_singularity "ON" \
    --p_singularity_bindpath ${singularity_bindpath} \
    --p_singularity_img ${singularity_img}
)
xgboost_jobid=$( echo $xgboost_jobid | awk '{print $4}' )

# GO-directness
${net_eval_dir}go_directness_sequential_corrected_pval \
    --p_in_dir_go "${evaluation_dir}/go250_xgboost_fixed_penalty_corrected_pval/" \
    --l_nbr_edges_per_reg 3 50 100 250 \
    --p_binding_event ${binding_label_stableid} \
    --p_out_eval "${evaluation_dir}/go250_xgboost_fixed_penalty_corrected_pval/eval_go250_directness_xgboost.tsv" \
    --p_out_logs "${log_dir}" \
    --data "go250_direct_${tissue}_xgboost" \
    --flag_slurm "ON" \
    --flag_singularity "ON" \
    --p_singularity_bindpath ${singularity_bindpath} \
    --p_singularity_img ${singularity_img} \
    --jobid_dependency ${xgboost_jobid}

