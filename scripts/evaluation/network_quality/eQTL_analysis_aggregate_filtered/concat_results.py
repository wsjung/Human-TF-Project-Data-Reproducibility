import pandas as pd
import numpy as np
import itertools
import argparse
import os
import re
import sys

parser = argparse.ArgumentParser(description="Take subsets of master df results, corresponding to each human tf network map")
parser.add_argument("--master_path", required=True, help="path to master networks, with variant support and count columns")
parser.add_argument("--variant_gene_pairs_path", required=True, help="path to GTEx eQTL files containing significant variant-gene pairs")
parser.add_argument("--variant_file_number", required=True, help="eqtl variant file number")
parser.add_argument("--variant_type", required=True, help="eqtl variant type. Should be '_filtered.tsv', '_unique.tsv', or '_ubiquitous.tsv'. Should also work with original data if you use '.v8.signif_variant_gene_pairs.txt'. Must put underline/period at beginning of string and file type at end of string")
parser.add_argument("--output_path", required=True, help="path to output human tf network maps with extra eQTL support and count columns")
args = parser.parse_args()

# Load variant file
files = sorted(file for file in os.listdir(args.variant_gene_pairs_path) if file.endswith(args.variant_type))
print(int(args.variant_file_number))
variant_file_name = files[int(args.variant_file_number)]
eqtl_tissue = os.path.basename(variant_file_name).partition(args.variant_type)[0]
eqtl_tissue = eqtl_tissue.lower()
print(f'eqtl tissue: {eqtl_tissue}')

df_master_final = pd.DataFrame()

# Add results across all batch df's to get final df_master
i = 1
for master_file in sorted(os.listdir(os.path.join(args.master_path, eqtl_tissue))):
    print(f'batch {i}')
    df_master = pd.read_csv(os.path.join(args.master_path, eqtl_tissue, master_file), sep="\t")
    # create indices based on tf and gene id's
    df_master = df_master.set_index(['TF','GENE','pred_prob'])
    if df_master_final.empty:
        df_master_final = df_master
    else:
        df_master_final = df_master_final + df_master
    i = i + 1

# reset the index to get rid of MultiIndex
df_master_final = df_master_final.reset_index()

# Save file in the directory
output_path = os.path.join(args.output_path, "final_output")
os.makedirs(output_path, exist_ok=True)
df_master_final.to_csv(os.path.join(output_path, f"{eqtl_tissue}.tsv"), sep="\t", index=False)

print('done')
