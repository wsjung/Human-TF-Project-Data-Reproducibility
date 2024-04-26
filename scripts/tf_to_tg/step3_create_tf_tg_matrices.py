import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="create TF_TG score matrices for each FANTOM5 sample")
parser.add_argument("--input_scores_path", required=True, help="path to annotated TF to promoter files for each sample, which contain tf_tg score columns")
parser.add_argument("--index", required=True, help="index corresponding to which 10 samples we're looking at")
parser.add_argument("--output_path", required=True, help="path to output tf_tg score matrices")
args = parser.parse_args()

scores_df_first_sample = pd.read_csv(os.path.join(args.input_scores_path, "0007.add_tf_tg_scores.txt"), sep="\t", comment="#")
unique_tf_ensg_ids = scores_df_first_sample.tf_ensg_id.unique()
unique_tg_ensg_ids = scores_df_first_sample.tg_ensg_id.unique()
zeros_matrix = np.zeros((len(unique_tf_ensg_ids),len(unique_tg_ensg_ids)))

i =  int(args.index)*10 + 7
stop = int(args.index)*10 + 17
stop_final = 1836

def add_to_matrix(row):
    tf_tg_score = row["tf_tg_score"]
    tf_ensg_id = row["tf_ensg_id"]
    tg_ensg_id = row["tg_ensg_id"]
    if(tf_tg_score > matrix_df.loc[tf_ensg_id, tg_ensg_id]):
        matrix_df.loc[tf_ensg_id, tg_ensg_id] = tf_tg_score

while((i < stop) and i < (stop_final)):
    scores_df = pd.read_csv(os.path.join(args.input_scores_path, f"{i:04d}.add_tf_tg_scores.txt"), sep="\t", comment="#")
    matrix_df = pd.DataFrame(zeros_matrix, columns = unique_tg_ensg_ids, index = unique_tf_ensg_ids)
    scores_df.apply(add_to_matrix, axis=1)
    matrix_df.to_csv(os.path.join(args.output_path, f"{i:04d}.create_tf_tg_matrices.txt"), sep="\t", index=True)
    print(i)
    i = i + 1

print("done")

