import pandas as pd
import numpy as np
import itertools
import argparse
import os
import re
import sys

parser = argparse.ArgumentParser(description="Determine which remap tfs are in the networks, so that we can create batches")
parser.add_argument("--master_path", required=True, help="path to master df")
parser.add_argument("--remap_path", required=True, help="path to remap files for each tf, containing tf, tg, and subregion of enhancer/promoter region associated with tg")
parser.add_argument("--output_path", required=True, help="path to tfs")
args = parser.parse_args()

# load master network map
df_master = pd.read_csv(args.master_path, sep="\t", header=None, names=['TF','GENE','pred_prob'])

# create file for output
f = open(os.path.join(args.output_path, 'all_tfs.txt'), "w")

# Loop through each remap file
for remap_file_name in sorted(os.listdir(args.remap_path)):
    # each remap file corresponds to a tf. See if that tf is in the master df at all
    tf = remap_file_name.partition(".bed")[0]
    df_master_filtered = df_master[df_master['TF']==tf]
    matching_tfs = df_master_filtered.shape[0]
    # if tf is in master df...
    if matching_tfs > 0:
        f.write(f'{remap_file_name}\n')

print('done')
