import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="add HGNC id column to TF to promoter file")
parser.add_argument("--tf_to_promoter_path", required=True, help="path to input Tf to promoter files")
parser.add_argument("--annotations_path", required=True, help="path to promoter-gene annotation file")
parser.add_argument("--file_name", required=True, help="TF to promoter file name, since the larger file was split into smaller files for efficiency")
parser.add_argument("--output_path", required=True, help="path to output TF to promoter file with HGNC id column")
args = parser.parse_args()

promoter_scores_path = os.path.join(args.tf_to_promoter_path, args.file_name)
promoter_scores_df = pd.read_csv(promoter_scores_path, sep="\t", comment="#",names=['chr_name','promoter_start','promoter_stop','promoter_id','promoter_score','promoter_strand','motif_id','tf_id','motif_start','motif_stop','motif_score','motif_strand'])

annotations_df = pd.read_csv(args.annotations_path, sep="\t", comment="#",names = ['CAGE_Peak_ID','Transcript_name','Distance','GeneID','HGNC/MGI_ID','UniProt_ID','Gene_name','Gene_symbol','Gene_synonyms','Gene_source'])

def get_gene_id(promoter_id):
    row = annotations_df.loc[annotations_df['CAGE_Peak_ID'] == promoter_id]
    if row.shape[0] > 1:
        raise Exception("only one HGNC ID should be associated with the promoter")
    hgnc_id = row.iloc[0]['HGNC/MGI_ID']
    print(hgnc_id)
    return(hgnc_id)

promoter_scores_df['hgnc_id'] = promoter_scores_df['promoter_id'].apply(get_gene_id)

promoter_scores_df.to_csv(os.path.join(args.output_path, args.file_name), sep="\t", index=False)

print("done")
