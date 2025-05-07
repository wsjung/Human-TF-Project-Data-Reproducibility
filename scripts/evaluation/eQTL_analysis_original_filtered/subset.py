import pandas as pd
import numpy as np
import itertools
import argparse
import os
import re
import sys

parser = argparse.ArgumentParser(description="Take subsets of master df results, corresponding to each human tf network map")
parser.add_argument("--human_network_files_path", required=True, help="path to human network files")
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

# load master network map, which has columns for eQTL support and count
df_master = pd.read_csv(os.path.join(args.master_path, eqtl_tissue, "master.tsv"), sep="\t")

for human_network_file in sorted(os.listdir(args.human_network_files_path)):
    df_human_network = pd.read_csv(os.path.join(args.human_network_files_path, human_network_file), sep="\t", header=None, names=['TF','GENE','pred_prob'])
    network_tissue = os.path.basename(human_network_file).partition("_xgboost")[0]
    print(f'network tissue: {network_tissue}')
    # make directory for network tissue if not made yet
    output_path = os.path.join(args.output_path, "final_output_per_human_network_tissue", network_tissue)
    os.makedirs(output_path, exist_ok=True)
    # filter master df based on edges in this specific tissue's human network map 
    df_master_filtered = df_master[pd.Series(list(zip(df_master['TF'], df_master['GENE']))).isin(list(zip(df_human_network['TF'], df_human_network['GENE'])))]
    # create indices based on tf and gene id's
    df_human_network = df_human_network.set_index(['TF','GENE'])
    df_master_filtered = df_master_filtered.set_index(['TF','GENE'])
    # add revelant columns that were determined from the above analysis (using the new index to align the rows properly)
    df_human_network['eqtl_support'] = df_master_filtered['eqtl_support']
    df_human_network['variant_count'] = df_master_filtered['variant_count']
    # reset the index to get rid of MultiIndex
    df_human_network = df_human_network.reset_index()
    # Save file in the directory
    df_human_network.to_csv(os.path.join(output_path, f"{eqtl_tissue}.tsv"), sep="\t", index=False)

print('done')
