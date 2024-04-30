import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="add tf_ensg_id column, hgnc_id column, entrezgene_id column, and tg_ensg_id column to TF to promoter file")
parser.add_argument("--tf_to_promoter_path", required=True, help="path to input TF to promoter file")
parser.add_argument("--activity_levels_path", required=True, help="path to promoter-gene activity levels file, with first 3 rows removed")
parser.add_argument("--fantom5_ensg_mapping_path", required=True, help="FANTOM5 gene ensg_id mapping file")
parser.add_argument("--lambert_tfs_ensg_mapping_path", required=True, help = "lambert tfs ensg_id mapping file")
parser.add_argument("--file_name", required=True, help="TF to promoter file name, since the larger file was split into smaller files for efficiency")
parser.add_argument("--output_path", required=True, help="path to output TF to promoter file with tf_ensg_id column, hgnc_id column, entrezgene_id column, and tg_ensg_id column")
parser.add_argument("--output_path_bad_tfs", required=True, help="path to tfs without ensg id mappings")
args = parser.parse_args()

promoter_scores_path = os.path.join(args.tf_to_promoter_path, args.file_name)
promoter_scores_df = pd.read_csv(promoter_scores_path, sep="\t", comment="#", names=['chr_name','promoter_start','promoter_stop','promoter_id','promoter_score','promoter_strand','motif_id','tf_id','motif_start','motif_stop','motif_score','motif_strand'])

activity_levels_df = pd.read_csv(args.activity_levels_path, sep="\t", comment="#", header=None)
# Create a dict using promoter_id, entrezgene_id, hgnc_id
temp_activity_levels_df = activity_levels_df.iloc[:,[0,4,5]]
# Let promoter_id be the key
temp_activity_levels_dict = temp_activity_levels_df.set_index(0).T.to_dict('list')

fantom5_ensg_mapping_df = pd.read_csv(args.fantom5_ensg_mapping_path, sep="\t", comment="#")
lambert_tfs_ensg_mapping_df = pd.read_csv(args.lambert_tfs_ensg_mapping_path, sep="\t", comment="#")

new_promoter_scores = [] # Used to create new promoter scores df with ensg mappings
tfs_without_ensg_id = []

def get_ensg_ids(promoter_scores_row):
    # Add tf_ensg_id column based on mapping for tf id
    promoter_scores_row['tf_id_upper'] = promoter_scores_row['tf_id'].upper()
    tf_id_upper = promoter_scores_row['tf_id_upper']
    tf_ensg_id_row = lambert_tfs_ensg_mapping_df.loc[lambert_tfs_ensg_mapping_df['gene_name'] == tf_id_upper]
    if tf_ensg_id_row.shape[0] == 0:
        tfs_without_ensg_id.append(tf_id_upper)
        return
    if tf_ensg_id_row.shape[0] > 1:
        raise Exception("only one ENSG ID should be associated with the tf")
    promoter_scores_row['tf_ensg_id'] = tf_ensg_id_row.iloc[0]['gene_id']
    # Add tg_ensg_id column based on mapping for tg id
    promoter_id = promoter_scores_row['promoter_id']
    hgnc_id = temp_activity_levels_dict[promoter_id][1]
    entrezgene_id = temp_activity_levels_dict[promoter_id][0]
    tg_ensg_id_count = 0
    if(not pd.isna(hgnc_id)): # If hgnc_id isn't na
        tg_ensg_ids = fantom5_ensg_mapping_df.loc[fantom5_ensg_mapping_df['gene_source'] == hgnc_id]
        tg_ensg_id_count = tg_ensg_ids.shape[0]
    elif(not pd.isna(entrezgene_id)): # If entrezgene_id isn't na
        tg_ensg_ids = fantom5_ensg_mapping_df.loc[fantom5_ensg_mapping_df['gene_source'] == entrezgene_id]
        tg_ensg_id_count = tg_ensg_ids.shape[0]
    for i in range(0, tg_ensg_id_count):
        temp_row = promoter_scores_row
        temp_row['entrezgene_id'] = entrezgene_id
        temp_row['hgnc_id'] = hgnc_id
        temp_row['tg_ensg_id'] = tg_ensg_ids.iloc[i]['gene_id']
        new_promoter_scores.append(temp_row.to_dict())

promoter_scores_df.apply(get_ensg_ids, axis=1) # Look at each row of promoter_scores
new_promoter_scores_df = pd.DataFrame(new_promoter_scores)
new_promoter_scores_df = new_promoter_scores_df.drop_duplicates()
new_promoter_scores_df.to_csv(os.path.join(args.output_path, args.file_name), sep="\t", index=False)

tfs_without_ensg_id_df = pd.DataFrame(tfs_without_ensg_id)
tfs_without_ensg_id_df = tfs_without_ensg_id_df.drop_duplicates()
tfs_without_ensg_id_df.to_csv(os.path.join(args.output_path_bad_tfs, args.file_name), sep="\t", index=False)
print("done")
