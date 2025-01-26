#!/bin/bash

# default value
tissue=""

usage() {
    echo "Usage: $0 <tissue>"
    exit 1
}

# parse args
while [[ "$#" -gt 0 ]]; do
    case $1 in
        *) tissue=$1 ;;
    esac
    shift
done
if [[ -z "$tissue" ]]; then
    usage
fi
echo "TISSUE ${tissue}"

input_dir=/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_v29_output/${tissue}/
output_dir=/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/genie3_2024_09_07/${tissue}/
mkdir -p ${output_dir}genie3_batches
head -q -n 1 ${input_dir}no_design_gene_counts_vst_wgcna_genie3.csv > ${output_dir}all_genes_comma_separated.csv
tr ',' '\n' < ${output_dir}all_genes_comma_separated.csv > ${output_dir}all_genes_newline_separated.txt
split -n l/199 -a 3 --numeric-suffixes=1 --additional-suffix=.txt ${output_dir}all_genes_newline_separated.txt ${output_dir}genie3_batches/batch

echo "genie3 batches saved to ${output_dir}genie3_batches"
