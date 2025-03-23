#!/bin/bash

p_src_data=/scratch/mblab/d.p.ruskin/human/data/protein_coding_high_expr_filter/tpm_concat_output/

tail -n +2 ${p_src_data}np3_tpm_concat.txt | cut -f1 > ${p_src_data}all_genes.txt
mkdir -p ${p_src_data}bart_batches
split -n l/200 -a 3 --numeric-suffixes=1 --additional-suffix=.txt ${p_src_data}all_genes.txt ${p_src_data}bart_batches/batch

echo 'batches saved to:'
echo ${p_src_data}bart_batches
