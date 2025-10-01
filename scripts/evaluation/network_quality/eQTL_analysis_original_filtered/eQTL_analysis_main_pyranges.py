import pandas as pd
import numpy as np
import itertools
import argparse
import os
import re
import sys
import pyranges as pr

parser = argparse.ArgumentParser(description="Determine tissue-specific eQTL support and variant count for subset of edges in master human tf network map")
parser.add_argument("--master_path", required=True, help="master human tf network map which is the union of edges across tissue types")
parser.add_argument("--remap_batch", required=True, help="list of remap tf's in current batch")
parser.add_argument("--remap_path", required=True, help="path to remap files for each tf, containing tf, tg, and subregion of enhancer/promoter region associated with tg")
parser.add_argument("--variant_gene_pairs_path", required=True, help="path to GTEx eQTL files containing significant variant-gene pairs")
parser.add_argument("--variant_file_number", required=True, help="eqtl variant file number")
parser.add_argument("--variant_type", required=True, help="eqtl variant type. Should be '_filtered.tsv', '_unique.tsv', or '_ubiquitous.tsv'. Should also work with original data if you use '.v8.signif_variant_gene_pairs.txt'. Must put underline/period at beginning of string and file type at end of string")
parser.add_argument("--gencode_ensg_mapping_path", required=True, help="gencode v29 gene ensg_id mapping file")
parser.add_argument("--output_path", required=True, help="path to output master human tf network map with extra eQTL support and count columns")
args = parser.parse_args()

# load master network map
df_master = pd.read_csv(args.master_path, sep="\t")
# add 2 columns to the master df
df_master['eqtl_support'] = 0
df_master['variant_count'] = 0

# Load variant file
files = sorted(file for file in os.listdir(args.variant_gene_pairs_path) if file.endswith(args.variant_type))
print(int(args.variant_file_number))
variant_file_name = files[int(args.variant_file_number)]
eqtl_tissue = os.path.basename(variant_file_name).partition(args.variant_type)[0]
eqtl_tissue = eqtl_tissue.lower()
print(f'eqtl tissue: {eqtl_tissue}')
df_variant_gene_pairs = pd.read_csv(os.path.join(args.variant_gene_pairs_path,variant_file_name), sep="\t")
df_variant_gene_pairs['Start'] = df_variant_gene_pairs['variant_id'].str.split('_', expand=True)[1].astype(int)
# Make 'End' column 1 + Start so that pyranges can count overlaps for the one 'Start' position (start is included and end is excluded in pyranges)
df_variant_gene_pairs['End'] = df_variant_gene_pairs['Start'] + 1
df_variant_gene_pairs['Chromosome'] = df_variant_gene_pairs['variant_id'].str.split('_', expand=True)[0]

# load mapping file and create dict
df_gencode_ensg_mapping = pd.read_csv(args.gencode_ensg_mapping_path, sep="\t")
df_gene_name_id_dict = {row['gene_name']:row['gene_id'] for _,row in df_gencode_ensg_mapping.iterrows()}

# load list of remap tfs in this batch
with open(args.remap_batch, 'r') as file:
    # read all the lines of the file into a list
    remap_files = file.read().splitlines()

# Loop through each remap subregion (subregions that TFs actually bind to within promoter/enhancers)
for remap_file_name in remap_files:
    print(remap_file_name)
    # each remap file corresponds to a tf. Look at where tf is in the master df
    tf = remap_file_name.partition(".bed")[0]
    df_master_filtered = df_master[df_master['TF']==tf]
    # load remap file and map target gene names to their ensg id's
    df_remap = pd.read_csv(os.path.join(args.remap_path, remap_file_name), sep="\t")
    df_remap['gene_ensg'] = df_remap['gene'].map(df_gene_name_id_dict)
    # edit Start and End columns to align with what pyranges needs
    df_remap['Start'] = df_remap['Start_overlap']
    df_remap['End'] = df_remap['End_overlap']
    # get list of the unique target genes for this tf. See how many of these tf-tg edges are in master df
    unique_genes = df_remap['gene_ensg'].unique()
    df_master_filtered = df_master_filtered[df_master_filtered['GENE'].isin(unique_genes)]
    # for each edge...
    for i in range(df_master_filtered.shape[0]):
        # filter remap file for just these edges
        gene = df_master_filtered.iloc[i]['GENE']
        df_remap_filtered = df_remap[df_remap['gene_ensg']==gene]
        # filter eQTL data to only include variants with a target gene of the edges
        df_variant_gene_pairs_filtered = df_variant_gene_pairs[df_variant_gene_pairs['gene_id']==gene]
        if df_variant_gene_pairs_filtered.shape[0] > 0:
            # create pyranges objects
            pr_remap_filtered = pr.PyRanges(df_remap_filtered)
            pr_variant_gene_pairs_filtered = pr.PyRanges(df_variant_gene_pairs_filtered)
            # count the exact number of variants inside promoter/enhancer regions
            pr_remap_filtered = pr_remap_filtered.count_overlaps(pr_variant_gene_pairs_filtered)
            count = pr_remap_filtered.NumberOverlaps.sum()
            if count > 0:
                df_master.loc[((df_master['TF'] == tf) & (df_master['GENE'] == gene)), 'eqtl_support'] = 1
                df_master.loc[((df_master['TF'] == tf) & (df_master['GENE'] == gene)), 'variant_count'] = count

# make directory for this eqtl tissue's batch output
output_path = os.path.join(args.output_path, "batch_outputs_per_eQTL_tissue", eqtl_tissue)
os.makedirs(output_path, exist_ok=True)
#df_master.to_csv(os.path.join(output_path, f"{os.path.basename(args.remap_batch).partition(".txt")[0]}.tsv"), sep="\t", index=False)
df_master.to_csv(os.path.join(output_path, os.path.basename(args.remap_batch)), sep="\t", index=False)

print('done')
