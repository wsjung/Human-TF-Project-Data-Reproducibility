import argparse
import os
import shutil
import pandas as pd

parser = argparse.ArgumentParser(description="renames tf_tg matrices and puts them in another output path, replacing sample numbers with the actual sample names")
parser.add_argument("--activity_levels_first_row_path", required=True, help="path to first row of activity levels file")
parser.add_argument("--output_with_sample_numbers_path", required=True, help="original output tf_tg score matrices, with sample numbers in the names")
parser.add_argument("--output_path", required=True, help="path to output tf_tg score matrices, this time with sample names instead of sample numbers")
args = parser.parse_args()

activity_levels_df = pd.read_csv(args.activity_levels_first_row_path, sep="\t", comment="#", header=None)

i = 7
stop = activity_levels_df.shape[1]
while(i < stop):
    shutil.copyfile(os.path.join(args.output_with_sample_numbers_path, f"{i:04d}.create_tf_tg_matrices.txt"), os.path.join(args.output_path, f"{activity_levels_df.iloc[0,i]}"))
    print(i)
    i = i + 1

print("done")
