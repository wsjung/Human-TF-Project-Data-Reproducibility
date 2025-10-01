import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="creates batches for LASSO")
parser.add_argument("--input_path", required=True, help="input path")
parser.add_argument("--tissue", required=True, help="tissue")
parser.add_argument("--aggregate", action='store_true', help="use if doing aggregate")
args = parser.parse_args()

print(args.tissue)

# (0) set variables
tissue_path = os.path.join(args.input_path, args.tissue)
bart_batches_path = os.path.join(tissue_path, "bart_batches")
lasso_batches_path = os.path.join(tissue_path, "lasso_batches")
if not os.path.exists(lasso_batches_path):
    os.makedirs(lasso_batches_path)

# (1) load data
df_tpm = pd.read_csv(os.path.join(tissue_path, "np3_tpm.txt"), sep="\t", index_col=0)
df_tpm_tfs = pd.read_csv(os.path.join(tissue_path, "np3_tpm_tfs.txt"), sep="\t", index_col=0)

# (2) subset df_tpm for batch genes
batch_num_max = 101
if args.aggregate:
    batch_num_max = 201

for batch_num in range(1,batch_num_max):
    batch_num = str(batch_num).zfill(3) # e.g. 1 -> 001
    batch_filename = f"batch{batch_num}.txt"
    print(batch_filename)

    # (2.1) load batch genes
    batch_genes = np.loadtxt(os.path.join(bart_batches_path, batch_filename), dtype=str)

    # (2.2) subset for TPM of batch genes
    df_tpm_batch_tgs = df_tpm[df_tpm.index.isin(batch_genes)]

    # (2.3) TPM of batch genes and TFs
    df_tpm_batch = pd.concat([df_tpm_batch_tgs, df_tpm_tfs]).drop_duplicates()

    # (2.4) save
    df_tpm_batch.to_csv(os.path.join(lasso_batches_path, batch_filename), sep="\t", index=True)

print(f"check {lasso_batches_path}")
print("DONE")
