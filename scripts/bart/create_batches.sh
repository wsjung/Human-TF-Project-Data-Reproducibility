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

p_src_data=/scratch/mblab/d.p.ruskin/human/GTEx_RNAseq/DATA/DTO_INPUT/protein_coding_high_expr_filter/GTEx_processing_gene_counts_tpm_v29_output/${tissue}/

tail -n +2 ${p_src_data}np3_tpm.txt | cut -f1 > ${p_src_data}all_genes.txt
mkdir -p ${p_src_data}bart_batches
split -n l/100 -a 3 --numeric-suffixes=1 --additional-suffix=.txt ${p_src_data}all_genes.txt ${p_src_data}bart_batches/batch

echo 'batches saved to:'
echo ${p_src_data}bart_batches
