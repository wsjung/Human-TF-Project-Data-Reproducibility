import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="extracts features")
parser.add_argument("--data_path", required=True, help="path to input data file")
parser.add_argument("--output_path", required=True, help="path to output features")
args = parser.parse_args()


# load data
df = pd.read_csv(args.data_path, sep="\t")
df = df.rename(columns={'TF':'REGULATOR','GENE':'TARGET'})


# for each feature column, extract into its own file
cols = df.columns[2:]
for col in cols:
    df_feat = df[['REGULATOR','TARGET',col]]
    df_feat.to_csv(os.path.join(args.output_path, f"net_{col}.tsv"), sep="\t", index=False, header=False)

# write the binding label file for evaluation
df_label = df[['REGULATOR','TARGET','LABEL']]
df_label_pos = df_label[df_label['LABEL']==1]
df_label_pos = df_label_pos.drop(columns=['LABEL'])
df_label_pos.to_csv(os.path.join(args.output_path, "binding_labels.txt"), sep="\t", index=False, header=True)


print("DONE")
