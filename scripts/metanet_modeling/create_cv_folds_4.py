import pandas as pd
import numpy as np
import os
import argparse
import random
from sklearn.model_selection import StratifiedKFold

parser = argparse.ArgumentParser(description="creates CV folds")
parser.add_argument("--tissue", required=True, help="tissue")
parser.add_argument("--base_path", required=True, help="base path")
parser.add_argument("--model", required=True, help="model (e.g. excl_perturbation)")
parser.add_argument("--binding_threshold", type=int, default=10, help="binding threshold (default=10)")
parser.add_argument("--num_CV_folds", type=int, default=10, help="number of CV folds")
parser.add_argument("--seed", type=int, default=42, help="random seed")
args = parser.parse_args()

## set random seed
#random.seed(args.seed)

# create dir to store CV data
output_path = os.path.join(args.base_path, "results", args.tissue, f"{args.binding_threshold}.0_threshold_{args.model}", "CV_folds")
if not os.path.exists(output_path):
    os.makedirs(output_path)

# load data
input_path = os.path.join(args.base_path, "input_data", args.tissue, f"{args.binding_threshold}.0_threshold_{args.model}", f"{args.tissue}_{args.binding_threshold}.0_threshold_{args.model}_nlogrank_minmax.txt")
df = pd.read_csv(input_path, sep="\t")

# preprocess (abs value)
for col in df.columns[3:]:
    df[col] = df[col].abs()

# split to features and response dfs
df_Y = df['LABEL']
df_X = df.drop(columns=['LABEL'])

# split CV folds
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=args.seed)

for fold, (train_idx, test_idx) in enumerate(skf.split(df_X, df_Y)):
    print(f"Fold {fold}")

    df_fold_train = df.iloc[train_idx,]
    df_fold_test = df.iloc[test_idx,]

    # save CV folds
    df_fold_train.to_csv(os.path.join(output_path, f"fold{fold}_train_data.txt"), sep="\t", index=False)
    df_fold_test.to_csv(os.path.join(output_path, f"fold{fold}_test_data.txt"), sep="\t", index=False)
    print(f"CV folds saved to {output_path}")

print("DONE")
