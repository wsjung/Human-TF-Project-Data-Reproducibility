import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="sums all paths from TF to TG")
parser.add_argument("--input_file", help="input file")
parser.add_argument("--output_file", help="output file")
parser.add_argument("--gencode", help="gencode mapping file")
args = parser.parse_args()

# read file
df = pd.read_csv(args.input_file, sep="\t")
df_gencode = pd.read_csv(args.gencode, sep="\t")
gencode_dict = { row['gene_name']:row['gene_id'] for _,row in df_gencode.iterrows() }

# map genes to ensg
df['gene_ensg'] = df['gene'].map(gencode_dict)

# sum up all TF-TG scores
df = df.groupby(['TF_ensg','gene_ensg'], as_index=False)['score'].sum()

# write to file
df = df.rename(columns={'TF_ensg':'TF', 'gene_ensg':'gene'})
df.to_csv(args.output_file, sep="\t", index=False)

print("DONE")
