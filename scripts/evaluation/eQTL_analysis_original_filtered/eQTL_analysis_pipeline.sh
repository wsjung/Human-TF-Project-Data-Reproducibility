#! /bin/bash


# FILE PATHS (MUST UPDATE BEFORE RUNNING THIS SCRIPT)

conda_env_dir=/scratch/mblab/d.p.ruskin/conda/human_network

# Need to create this empty directory before running this script
pipeline_dir=/scratch/mblab/d.p.ruskin/human/data/xgboost_all_gtex_tissues/eQTL_analysis_pipeline_original_filtered/

scripts_dir=/scratch/mblab/d.p.ruskin/human/scripts/eQTL_analysis_original_filtered/

# Must manually copy the 36 human network files into one folder before running this script
human_network_files_dir=/scratch/mblab/d.p.ruskin/human/data/xgboost_all_gtex_tissues/10.0_threshold_lasso_excl_genie3_wgcna/

remap_dir=/scratch/mblab/d.p.ruskin/human/data/new_REMAP_binding_labels/Remap_Data_TFWise_ENSG_1_qval0.01_mapping_annotated_overlapping_coordinates/

variant_gene_pairs_dir=/scratch/mblab/jungw/human_TF_project/data/xgboost_all_gtex_tissues/eval_tissue_specificity/eqtl_filtered/

# Should be '_filtered.tsv', '_unique.tsv', or '_ubiquitous.tsv'. Also works with original data if you use '.v8.signif_variant_gene_pairs.txt', as long as the files are already unzipped. Must put underline/period at beginning of string and file type at end of string
variant_type="_filtered.tsv"

gencode_ensg_mappings=/scratch/mblab/d.p.ruskin/human/data/gencode_v29/gencode.v29.protein_coding_genes_aliases.tsv


# SCRIPTS (NO NEED TO UPDATE)

concat_network_maps_job_id=$(sbatch --parsable ${scripts_dir}concat_network_maps.sh \
	${conda_env_dir} \
	${pipeline_dir} \
	${scripts_dir} \
	${human_network_files_dir})
echo 'concat_network_maps_job_id'
echo $concat_network_maps_job_id


find_remap_tfs_job_id=$(sbatch --dependency=afterok:$concat_network_maps_job_id --parsable ${scripts_dir}find_remap_tfs.sh \
	${conda_env_dir} \
	${pipeline_dir} \
	${scripts_dir} \
	${remap_dir})
echo 'find_remap_tfs_job_id'
echo $find_remap_tfs_job_id


slurmids=""
for i in {0..35};
do
        slurmids="$slurmids:$(sbatch --dependency=afterok:$find_remap_tfs_job_id --parsable ${scripts_dir}eQTL_analysis_main_pyranges.sh \
		"$i" \
		${conda_env_dir} \
		${pipeline_dir} \
		${scripts_dir} \
		${remap_dir} \
		${variant_gene_pairs_dir} \
		${variant_type} \
		${gencode_ensg_mappings})"
done
echo 'slurmids'
echo $slurmids


concat_results_job_id=$(sbatch --dependency=afterok$slurmids --parsable ${scripts_dir}concat_results.sh \
	${conda_env_dir} \
	${pipeline_dir} \
	${scripts_dir} \
	${variant_gene_pairs_dir} \
	${variant_type})
echo 'concat_results_job_id'
echo $concat_results_job_id


subset_job_id=$(sbatch --dependency=afterok:$concat_results_job_id --parsable ${scripts_dir}subset.sh \
	${conda_env_dir} \
	${pipeline_dir} \
	${scripts_dir} \
	${human_network_files_dir} \
	${variant_gene_pairs_dir} \
	${variant_type})
echo 'subset_job_id'
echo $subset_job_id
