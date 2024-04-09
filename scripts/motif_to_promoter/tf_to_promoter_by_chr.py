import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="groups motif hits by promoter region")
parser.add_argument("--bed_path", required=True, help="path to bed files directory, sorted by chromosome")
parser.add_argument("--fimo_path", required=True, help="path to fimo with overlaps files directory, sorted by chromosome")
parser.add_argument("--chromosome_strand", required=True, help="chromosome number")
parser.add_argument("--output_path", required=True, help="path to output tf_to_promoter files, sorted by chromosome")
args = parser.parse_args()

bed_path = os.path.join(args.bed_path, args.chromosome_strand + '.bed')
bed_df = pd.read_csv(bed_path, sep="\t", comment='#', names = ['chr_name','start','stop','ID','a','strand','b','c','d']) # Don't care about columns a, b, c, d
fimo_path = os.path.join(args.fimo_path, args.chromosome_strand + '.txt')
fimo_df = pd.read_csv(fimo_path, sep="\t", comment='#', names = ['motif_id','motif_alt_id','strand','score','p-value','q-value','matched_sequence','chr_name','sequence_start','sequence_stop'])
rows = []
i = 0
while i < bed_df.shape[0]:
    j = 0
    while(j < fimo_df.shape[0]):
        # If bed region is completely within fimo region
        if fimo_df.iloc[j]['sequence_start'] <= bed_df.iloc[i]['start'] and fimo_df.iloc[j]['sequence_stop'] >= bed_df.iloc[i]['stop']:
            promoter_score = fimo_df.iloc[j]['score'] * (bed_df.iloc[i]['stop']-bed_df.iloc[i]['start'] + 1)
            rows.append([bed_df.iloc[i]['chr_name'], bed_df.iloc[i]['start'], bed_df.iloc[i]['stop'], bed_df.iloc[i]['ID'], promoter_score, bed_df.iloc[i]['strand'], fimo_df.iloc[j]['motif_id'], fimo_df.iloc[j]['motif_alt_id'], fimo_df.iloc[j]['sequence_start'], fimo_df.iloc[j]['sequence_stop'], fimo_df.iloc[j]['score'], fimo_df.iloc[j]['strand']])
        # If fimo region partially overlaps with left side of bed region
        elif fimo_df.iloc[j]['sequence_start'] <= bed_df.iloc[i]['start'] and fimo_df.iloc[j]['sequence_stop'] >= bed_df.iloc[i]['start']:
            promoter_score = fimo_df.iloc[j]['score'] * (fimo_df.iloc[j]['sequence_stop'] - bed_df.iloc[i]['start'] + 1) 
            rows.append([bed_df.iloc[i]['chr_name'], bed_df.iloc[i]['start'], bed_df.iloc[i]['stop'], bed_df.iloc[i]['ID'], promoter_score, bed_df.iloc[i]['strand'], fimo_df.iloc[j]['motif_id'], fimo_df.iloc[j]['motif_alt_id'], fimo_df.iloc[j]['sequence_start'], fimo_df.iloc[j]['sequence_stop'], fimo_df.iloc[j]['score'], fimo_df.iloc[j]['strand']]) 
        # If fimo region partially overlaps with right side of bed region
        elif fimo_df.iloc[j]['sequence_start'] <= bed_df.iloc[i]['stop'] and fimo_df.iloc[j]['sequence_stop'] >= bed_df.iloc[i]['stop']:
            promoter_score = fimo_df.iloc[j]['score'] * (bed_df.iloc[i]['stop'] - fimo_df.iloc[j]['sequence_start'] + 1)
            rows.append([bed_df.iloc[i]['chr_name'], bed_df.iloc[i]['start'], bed_df.iloc[i]['stop'], bed_df.iloc[i]['ID'], promoter_score, bed_df.iloc[i]['strand'], fimo_df.iloc[j]['motif_id'], fimo_df.iloc[j]['motif_alt_id'], fimo_df.iloc[j]['sequence_start'], fimo_df.iloc[j]['sequence_stop'], fimo_df.iloc[j]['score'], fimo_df.iloc[j]['strand']])
        # If fimo region is completely within bed region
        elif fimo_df.iloc[j]['sequence_start'] >= bed_df.iloc[i]['start'] and fimo_df.iloc[j]['sequence_stop'] <= bed_df.iloc[i]['stop']:
            promoter_score = fimo_df.iloc[j]['score'] * (fimo_df.iloc[j]['sequence_stop'] - fimo_df.iloc[j]['sequence_start'] + 1)
            rows.append([bed_df.iloc[i]['chr_name'], bed_df.iloc[i]['start'], bed_df.iloc[i]['stop'], bed_df.iloc[i]['ID'], promoter_score, bed_df.iloc[i]['strand'], fimo_df.iloc[j]['motif_id'], fimo_df.iloc[j]['motif_alt_id'], fimo_df.iloc[j]['sequence_start'], fimo_df.iloc[j]['sequence_stop'], fimo_df.iloc[j]['score'], fimo_df.iloc[j]['strand']])
        j = j + 1
    i = i + 1
if rows:
    output_df = pd.DataFrame(rows, columns=('chr_name','promoter_start','promoter_stop','promoter_id','promoter_score','promoter_strand','motif_id','tf_id','motif_start','motif_stop','motif_score','motif_strand'))
    output_df = output_df.drop_duplicates()
    output_df.to_csv(os.path.join(args.output_path, f"{args.chromosome_strand}.tf_to_promoter_by_chr.txt"), sep="\t", index=False)

print("done")
