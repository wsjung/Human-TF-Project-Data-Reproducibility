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

output_dir=/scratch/mblab/d.p.ruskin/human/data/NP3_results_marbach/tissue_specific_gtex_tpm/${tissue}/bart_output/

# Create a tab-separated file with all the gene names
head -n1 -q ${output_dir}/batches_output/* > ${output_dir}headers.txt
tr '\n' '\t' < ${output_dir}headers.txt > ${output_dir}tab_separated_headers.tsv

arr=(${output_dir}batches_output/*)

file="${arr[0]}"

# keep first column of first batch (tf names). Add all the other batches without their first column
for f in "${arr[@]:1}"; do
   paste "$file" <(cut -d$'\t' -f2- "$f") > ${output_dir}_file.tmp && mv ${output_dir}_file.tmp ${output_dir}file.tmp
   file=${output_dir}file.tmp
done

tail -n +2 ${output_dir}file.tmp > ${output_dir}bottom_half.tsv # remove header because it's lined up weirdly

awk '{print}' ${output_dir}tab_separated_headers.tsv ${output_dir}bottom_half.tsv > ${output_dir}new_file.tmp # prepend with header now

#awk -F'\t' 'BEGIN { OFS = FS }; NF { NF -= 1 }; 1' <new_file.tmp> ${output_dir}net_bart.tsv # remove last column for consistency with np3 bart

sed 's/[ \t]\+$//' ${output_dir}new_file.tmp > ${output_dir}net_bart.tsv

echo 'results saved to:'
echo ${output_dir}net_bart.tsv
