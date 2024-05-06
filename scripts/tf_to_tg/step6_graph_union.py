import pandas as pd
import numpy as np
import argparse
import os
import re

parser = argparse.ArgumentParser(description="group fantom5 sample matrices by gtex tissue, and take the graph union of each tissue")
parser.add_argument("--sample_matrices_path", required=True, help="path to input matrices for each fantom5 sample, mapped to a tissue. Files are named {tissue}_{sample_name}")
parser.add_argument("--tissue", required=True, type = str, help="tissue we're focusing on")
parser.add_argument("--output_path", required=True, help="path to output matrices for each tissue")
args = parser.parse_args()

# Grab the sample matrices corresponding to the tissue (their file name contains the tissue name as a prefix)
tissue_samples = [filename for filename in os.listdir(args.sample_matrices_path) if filename.startswith(f"{args.tissue}_")]
print(f'Number of {args.tissue} samples: {len(tissue_samples)}')

# Read and concat the samples
sample_dfs = pd.concat(pd.read_csv(os.path.join(args.sample_matrices_path, sample_name), index_col=0, sep="\t", comment="#") for sample_name in tissue_samples)

# Take graph union
union_df = sample_dfs.groupby(level=0).max()

# Replace all special characters including whitespaces with underscores
gtex_tissue_underscore = re.sub(r'[^\w.]', '_', args.tissue)

union_df.to_csv(os.path.join(args.output_path, f"{gtex_tissue_underscore}.txt"), sep="\t", index=True)

print("done")
