import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="takes max peak score for overlapping remap peaks and calculates sum of scores for promoter/enhancer")
parser.add_argument("--input_path", required=True, help="input file with mapped peaks to promoter/enhancer")
parser.add_argument("--output_path", required=True, help="output path")
args = parser.parse_args()


def sum_max_score_by_basepair(df):

    sum_score = 0

    # for each basepair
    for pos in range(df['Start'].iloc[0], df['End'].iloc[0]):
        # peaks covering basepair
        df_peaks = df[(df['Start_remap'] <= pos) & (df['End_remap'] >= pos)]
        # max peak score at this bp
        max_score_pos = max(0, df_peaks['-log10q'].max())
        # sum
        sum_score = sum_score + max_score_pos

    return sum_score


# read files
print("reading files")
df = pd.read_csv(args.input_path, sep="\t")

if len(df) < 1:
    print("empty df")
    exit(1)

# preprocess
print("preprocessing")
df['ele_ID'] = df['Chromosome'].astype(str) + ":" + df['Start'].astype(str) + "-" + df['End'].astype(str)
tf = df['TF'].iloc[0]


# taking max score for each basepair
print("taking max score")
print(f"{df['ele_ID'].nunique()} unique promoter/enhancer regions")

list_ele_id = df['ele_ID'].unique()
list_scores = []
for ele_id in list_ele_id:
    df_ele = df[df['ele_ID'] == ele_id]
    score = sum_max_score_by_basepair(df_ele)
    list_scores.append(score)

# postprocess
print("postprocessing")
df_out = pd.DataFrame({'ele_ID':list_ele_id, 'score':list_scores})
df_out['Chromosome'] = df_out['ele_ID'].str.split(":", expand=True)[0]
df_out['Start'] = df_out['ele_ID'].str.split(":", expand=True)[1].str.split("-", expand=True)[0]
df_out['End'] = df_out['ele_ID'].str.split(":", expand=True)[1].str.split("-", expand=True)[1]
df_out['TF'] = tf
df_out = df_out.drop(columns=['ele_ID'])
df_out = df_out[['Chromosome','Start','End','TF','score']]

# write
print("writing")
df_out.to_csv(args.output_path, sep="\t", index=False)


print("DONE")
