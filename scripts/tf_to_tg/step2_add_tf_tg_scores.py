import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="add promoter_tg and tf_tg score columns to annotated TF to promoter file, creating a separate file for each FANTOM5 sample")
parser.add_argument("--tf_to_promoter_path", required=True, help="path to input annotated TF to promoter file")
parser.add_argument("--activity_levels_path", required=True, help="path to promoter-gene activity levels file, with first 3 rows removed")
parser.add_argument("--index", required=True, help="index corresponding to which 10 samples we're looking at")
parser.add_argument("--output_path", required=True, help="path to output annotated TF to promoter files with new promoter_tg and tf_tg score columns")
args = parser.parse_args()

#annotated_promoter_scores_path = os.path.join(args.tf_to_promoter_path, args.file_name)
annotated_promoter_scores_df = pd.read_csv(args.tf_to_promoter_path, sep="\t", comment="#")

activity_levels_df = pd.read_csv(args.activity_levels_path, sep="\t", comment="#", header=None)
def get_promoter_tg_score(promoter_id):
    #row = activity_levels_df[activity_levels_df[0] == promoter_id]
    #if row.shape[0] == 0:
    #    raise Exception("only one row should be associated with the promoter")
    #promoter_tg_score = row.iloc[0][i]
    promoter_tg_score = temp_activity_levels_dict[promoter_id][0]
    return(promoter_tg_score)

i =  int(args.index)*10 + 7
stop = int(args.index)*10 + 17
stop_final = activity_levels_df.shape[1]
while((i < stop) and i < (stop_final)):
    temp_activity_levels_df = activity_levels_df.iloc[:,[0,i]]
    temp_activity_levels_dict = temp_activity_levels_df.set_index(0).T.to_dict('list')
    temp_annotated_promoter_scores_df = annotated_promoter_scores_df
    temp_annotated_promoter_scores_df['promoter_tg_score'] = np.vectorize(get_promoter_tg_score)(temp_annotated_promoter_scores_df['promoter_id'])
    #temp_annotated_promoter_scores_df['promoter_id'].apply(get_promoter_tg_score)
    temp_annotated_promoter_scores_df['tf_tg_score'] = temp_annotated_promoter_scores_df['promoter_score']*temp_annotated_promoter_scores_df['promoter_tg_score']
    temp_annotated_promoter_scores_df.to_csv(os.path.join(args.output_path, f"{i:04d}.add_tf_tg_scores.txt"), sep="\t", index=False)
    i = i + 1
    print(i)

print("done")
