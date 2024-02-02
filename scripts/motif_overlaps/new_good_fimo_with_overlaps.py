import pandas as pd
import numpy as np
import argparse
import os

def get_chr_name(sequence):
    char_num1_num2 = sequence.replace(':',' ').replace('-',' ').split()
    return char_num1_num2[0]

def get_sequence_start(sequence):
    char_num1_num2 = sequence.replace(':',' ').replace('-',' ').split()
    return int(char_num1_num2[1])

parser = argparse.ArgumentParser(description="calculates motif hit scores based on overlaps")
parser.add_argument("--fimo_path", required=True, help="path to fimo output file")
parser.add_argument("--tf_ID", required=True, help="TF MEME ID")
parser.add_argument("--tf", required=True, help="TF")
parser.add_argument("--output_path", required=True, help="path to output directory")
args = parser.parse_args()


fimo_path = os.path.join(args.fimo_path, args.tf_ID, "fimo.tsv")

# for each lambert TF, check if overlaps exist and calculate corresponding scores

overlap = []  # List for keeping track of consecutive overlaps
new_fimo = []  # List for storing new fimo output, accounting for overlaps

try:
    df = pd.read_csv(fimo_path, sep="\t", comment='#')
    df['chr_name'] = df['sequence_name'].apply(get_chr_name)
    df['sequence_start'] = df['sequence_name'].apply(get_sequence_start) + df['start'] - 1
    df['sequence_stop'] = df['sequence_name'].apply(get_sequence_start) + df['stop'] - 1
    df = df.drop(['sequence_name','start','stop'], axis=1)
    df = df.sort_values(['chr_name', 'strand','sequence_start'])

    i = 1
    while i < df.shape[0]:
        if df.iloc[i]['chr_name'] == df.iloc[i-1]['chr_name'] and df.iloc[i]['strand'] == df.iloc[i-1]['strand'] and df.iloc[i]['sequence_start'] <= df.iloc[i-1]['sequence_stop']:
            score_A = df.iloc[i-1]['score']
            score_B = df.iloc[i]['score']
            if not overlap: # If overlap list is empty (previous motif hit was not overlapping the one before it)
                if score_A >= score_B:
                    #df.iloc[i,df.columns.get_loc('score')] = round((score_A * (df.iloc[i-1]['stop'] - df.iloc[i-1]['start']) + score_B * (df.iloc[i]['stop'] - df.iloc[i-1]['stop']))/(df.iloc[i]['stop'] - df.iloc[i-1]['start']), 4)
                    overlap.extend([score_A] * (df.iloc[i-1]['sequence_stop'] - df.iloc[i-1]['sequence_start']))
                    overlap.extend([score_B] * (df.iloc[i]['sequence_stop'] - df.iloc[i-1]['sequence_stop']))
                else:
                    #df.iloc[i,df.columns.get_loc('score')] = round((score_A * (df.iloc[i]['start'] - df.iloc[i-1]['start']) + score_B * (df.iloc[i]['stop'] - df.iloc[i]['start']))/(df.iloc[i]['stop'] - df.iloc[i-1]['start']),4)
                    overlap.extend([score_A] * (df.iloc[i]['sequence_start'] - df.iloc[i-1]['sequence_start']))
                    overlap.extend([score_B] * (df.iloc[i]['sequence_stop'] - df.iloc[i]['sequence_start']))
            else: # If overlap list is not empty (previous motif hit was overlapping the one before it)
                for x in range(len(overlap) - (df.iloc[i-1]['sequence_stop'] - df.iloc[i]['sequence_start']), len(overlap)):
                    if score_B > overlap[x]:
                        overlap[x] = score_B
                overlap.extend([score_B] * (df.iloc[i]['sequence_stop'] - df.iloc[i-1]['sequence_stop']))
            if df.iloc[i-1]['p-value'] < df.iloc[i]['p-value']:
                df.iloc[i,df.columns.get_loc('p-value')] = df.iloc[i-1]['p-value']
            if df.iloc[i-1]['q-value'] < df.iloc[i]['q-value']:
                df.iloc[i,df.columns.get_loc('q-value')] = df.iloc[i-1]['q-value']
            matched_sequence_a = df.iloc[i-1]['matched_sequence']
            matched_sequence_b = df.iloc[i]['matched_sequence']
            concatenate_from_here = df.iloc[i-1]['sequence_stop'] - df.iloc[i]['sequence_start'] + 1
            if df.iloc[i]['strand'] == "+":
                df.iloc[i,df.columns.get_loc('matched_sequence')] = matched_sequence_a + matched_sequence_b[concatenate_from_here:]
            elif df.iloc[i]['strand'] == "-":
                df.iloc[i,df.columns.get_loc('matched_sequence')] = matched_sequence_b + matched_sequence_a[concatenate_from_here:]
            df.iloc[i,df.columns.get_loc('sequence_start')] = df.iloc[i-1]['sequence_start']
            i = i + 1
        else:
            if overlap: # If there had been overlaps, then calculate total overlap score
                df.iloc[i-1,df.columns.get_loc('score')] = round(sum(overlap)/len(overlap),4) # Score is average of list
                overlap = [] # Set overlap list back to empty after calculating overlap score
            new_fimo.append(df.iloc[i-1].to_dict()) # Add motif hit to new fimo output
            i = i + 1
    # If the end of fimo has been reached, then add the final 
    if overlap: # If there had been overlaps, then calculate total overlap score
        df.iloc[i-1,df.columns.get_loc('score')] = round(sum(overlap)/len(overlap),4) # Score is average of list
        overlap = [] # Set overlap list back to empty after calculating overlap score
    new_fimo.append(df.iloc[i-1].to_dict()) # Add motif hit to new fimo output
except pd.errors.EmptyDataError:
    print("> error (empty data): no results!")

if new_fimo:
    new_fimo_df = pd.DataFrame(new_fimo)
    new_fimo_df = new_fimo_df.drop_duplicates()
    new_fimo_df.to_csv(os.path.join(args.output_path, f"{args.tf_ID}.new_good_fimo_with_overlaps.txt"), sep="\t", index=False)

print("done")
