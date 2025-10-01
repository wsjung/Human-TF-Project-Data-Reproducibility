import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="filters out non-significant (qval > 0.01) peaks")
parser.add_argument("--input_path", required=True, help="input path")
parser.add_argument("--filename", required=True, help="filename")
parser.add_argument("--output_path", required=True, help="output path")
args = parser.parse_args()

print(args.filename)

# read file
print("reading file")
df = pd.read_csv(os.path.join(args.input_path, args.filename), sep="\t")

# filter for significant peaks
print("filtering")
df = df[df['-log10q'] >= 2]

# rename cols
print("renaming cols")
df = df.rename(columns={'chr':'Chromosome', 'start':'Start', 'end':'End'})

# save
print("saving file")
df.to_csv(os.path.join(args.output_path, args.filename), sep="\t", index=False)

print("DONE")
