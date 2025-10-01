import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="annotates max scored promoter/enhancer regions with genes")
parser.add_argument("--input_path", required=True, help="input directory")
parser.add_argument("--filename", required=True, help='filename')
parser.add_argument("--TF_ensg", required=True, help="TF ENSG ID")
parser.add_argument("--gene_annotation", required=True, help="fantom5 gene annotation")
parser.add_argument("--output_path", required=True, help="output directory")
args = parser.parse_args()


# read in file
print("reading file")
try:
    df = pd.read_csv(os.path.join(args.input_path, args.filename), sep="\t")
except:
    print("file missing")
    exit(1)
df_annotation = pd.read_csv(args.gene_annotation, sep="\t")

df_merged = df.merge(df_annotation, on=['Chromosome','Start','End'], how='inner')
df_merged['TF_ensg'] = args.TF_ensg

# write
print("writing file")
df_merged.to_csv(os.path.join(args.output_path, args.filename), sep="\t", index=False)

print("DONE")
