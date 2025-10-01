import pandas as pd
import numpy as np
import os
import argparse
import pyranges as pr

parser = argparse.ArgumentParser(description="maps remap peaks to fantom5 elements")
parser.add_argument("--fantom5_path", required=True, help="path to fantom5 elements")
parser.add_argument("--remap_path", required=True, help="path to remap peak")
parser.add_argument("--output_path", required=True, help="output path")
args = parser.parse_args()


# read files
print("reading files")
df_fantom5 = pd.read_csv(args.fantom5_path, sep="\t")
df_remap = pd.read_csv(args.remap_path, sep="\t")

# preprocess
print("preprocessing")
## (i) col names Chromosome, Start, End
df_fantom5 = df_fantom5.rename(columns={'chr':'Chromosome','start':'Start','end':'End'})
df_remap = df_remap.rename(columns={'chr':'Chromosome','start':'Start','end':'End'})

## (ii) Start < End
df_fantom5['Start'], df_fantom5['End'] = df_fantom5[['Start', 'End']].min(axis=1), df_fantom5[['Start', 'End']].max(axis=1)
df_remap['Start'], df_remap['End'] = df_remap[['Start', 'End']].min(axis=1), df_remap[['Start', 'End']].max(axis=1)


# intersect
print("intersecting")
pr_fantom5 = pr.PyRanges(df_fantom5)
pr_remap = pr.PyRanges(df_remap)
pr_overlap = pr_fantom5.join(pr_remap, report_overlap=True, suffix="_remap")

# only keep overlaps by >=50% of peak region
print("keep only >=50% peak region overlaps")
df_overlap = pr_overlap.df
df_overlap['remap_length'] = df_overlap['End_remap'] - df_overlap['Start_remap']
df_overlap = df_overlap[['Chromosome','Start','End','Start_remap','End_remap','Overlap','remap_length','-log10q','TF']]
df_overlap = df_overlap[df_overlap['Overlap'] >= 0.5*df_overlap['remap_length']]

# save
print("writing")
df_overlap.to_csv(args.output_path, sep="\t", index=False)


print("DONE")
