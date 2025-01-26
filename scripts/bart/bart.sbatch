#!/bin/bash

#SBATCH --array=1-100
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/scratch/mblab/d.p.ruskin/human/logs/bart_%A_%a.out
#SBATCH -e /scratch/mblab/d.p.ruskin/human/logs/bart_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.p.ruskin@wustl.edu

#MUST SET THIS STRING BEFORE EACH TISSUE RUN
tissue="whole_blood"

# Data paths
p_src_data=/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_tpm_v29_output/${tissue}/
batch_dir=${p_src_data}/bart_batches
# Get batch file corresponding to the current task
batch_file=$(ls $batch_dir | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Code paths
p_src_code=/scratch/mblab/d.p.ruskin/NetProphet_3.0/
p_singularity_img=${p_src_code}singularity/s_np3.sif
p_scripts=/scratch/mblab/d.p.ruskin/human/scripts/bart/

# Output paths
p_out_dir=/scratch/mblab/d.p.ruskin/human/data/NP3_results_marbach/tissue_specific_gtex_tpm/${tissue}/bart_output/batches_output
mkdir -p ${p_out_dir}

# Load R environment
source ${p_src_code}src/helper/load_singularity.sh
export SINGULARITY_BINDPATH=/scratch/mblab

# Run the R script
singularity exec ${p_singularity_img} \
	Rscript --no-save --vanilla ${p_scripts}/better_build_net_bart.R \
	--p_in_expr_target ${p_src_data}/np3_tpm.txt \
     	--p_in_expr_reg ${p_src_data}/np3_tpm_tfs.txt \
     	--ntree 50 \
     	--p_out_dir ${p_out_dir} \
     	--p_src_code ${p_src_code} \
	--p_scripts ${p_scripts} \
	--batch_file $batch_dir/$batch_file

