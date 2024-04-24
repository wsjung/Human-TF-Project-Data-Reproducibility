import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="sums promoter scores for multiple motif hits in one promoter region, one chromosome strand at a time")
parser.add_argument("--promoter_scores_path", required=True, help="path to tf_to_promoter_by_chr output files")
parser.add_argument("--chromosome_strand", required=True, help="chr number")
parser.add_argument("--output_path", required=True, help="path to output directory")
args = parser.parse_args()

promoter_scores_path = os.path.join(args.promoter_scores_path, f"{args.chromosome_strand}.tf_to_promoter_by_chr.txt")

sum_score = 0.0 # Number for summing up promoter scores of motif hits in a promoter region
new_promoter_scores = []  # List for storing new tf_to_promoter scores, which sum up promoter scores for each motif hit in each promoter region

try:
    df = pd.read_csv(promoter_scores_path, sep="\t", comment='#')
    df = df.sort_values(['promoter_id', 'tf_id'])
    i = 1
    while i < df.shape[0]:
        if df.iloc[i]['promoter_id'] == df.iloc[i-1]['promoter_id'] and df.iloc[i]['tf_id'] == df.iloc[i-1]['tf_id']:
            if sum_score == 0:
                sum_score = sum_score + df.iloc[i-1]['promoter_score']
                sum_score = sum_score + df.iloc[i]['promoter_score']
            else:
                sum_score = sum_score + df.iloc[i]['promoter_score']
        else:
            if sum_score != 0: # If there were multiple hits in a region, then put the summed score in the df
                df.iloc[i-1,df.columns.get_loc('promoter_score')] = sum_score
                sum_score = 0.0 # Set sum_score back to 0 after it's not needed for this motif and promoter pair anymore
            new_promoter_scores.append(df.iloc[i-1].to_dict()) # Add row to new_promoter_scores output
        i = i + 1
    # If the end of fimo has been reached, then add the final 
    if sum_score != 0:
        df.iloc[i-1,df.columns.get_loc('promoter_score')] = sum_score
        sum_score = 0.0
    new_promoter_scores.append(df.iloc[i-1].to_dict())
except pd.errors.EmptyDataError:
    print("> error (empty data): no results!")

if new_promoter_scores:
    new_promoter_scores_df = pd.DataFrame(new_promoter_scores)
    new_promoter_scores_df = new_promoter_scores_df.drop_duplicates()
    new_promoter_scores_df.to_csv(os.path.join(args.output_path, f"{args.chromosome_strand}.sum_hits_by_chr.txt"), sep="\t", index=False)

print("done")
