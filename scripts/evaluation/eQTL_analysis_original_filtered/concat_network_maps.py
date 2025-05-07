import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="create a master network map which is the union of edges across tissue types")
parser.add_argument("--human_network_files_path", required=True, help="path to human network files")
parser.add_argument("--output_path", required=True, help="path to output")
args = parser.parse_args()

# create emtpty df, to store master network map
df_master = pd.DataFrame()
# load network maps and create a master network map which is the union of edges across tissue types
for human_network_file in sorted(os.listdir(args.human_network_files_path)):
    # only load tf and gene columns (df has no header so create TF and GENE column headers)
    df_human_network_edges = pd.read_csv(os.path.join(args.human_network_files_path, human_network_file), sep="\t", header=None, usecols=[0,1], names=['TF','GENE'])
    network_tissue = os.path.basename(human_network_file).partition("_xgboost")[0]
    print(f'network tissue: {network_tissue}')
    df_master = pd.concat([df_master, df_human_network_edges])
# drop duplicates, to get rid of duplicate edges (use ignore_index=True to avoid duplicate indices)
df_master = df_master.drop_duplicates(ignore_index=True)

df_master.to_csv(os.path.join(args.output_path, "xgboost_master.tsv"), sep="\t", index=False)

print('done')
