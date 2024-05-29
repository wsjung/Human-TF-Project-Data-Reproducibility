import argparse
import os
import shutil
import pandas as pd

parser = argparse.ArgumentParser(description="Adds GTEx tissue mapping to the start of each tf_tg matrix file name. Each matrix will now be named '{tissue}_{sample}'")
parser.add_argument("--fantom5_tissue_mapping_path", required=True, help="FANTOM5 sample to tissue mapping file")
parser.add_argument("--output_with_sample_names_path", required=True, help="tf_tg score matrices, named according to their FANTOM5 sample")
parser.add_argument("--output_path", required=True, help="path to output tf_tg score matrices, each named '{tissue}_{sample}'")
args = parser.parse_args()

fantom5_tissue_mapping_df = pd.read_csv(args.fantom5_tissue_mapping_path, sep="\t", comment="#")
# Create a dict using sample name and tissue
temp_fantom5_tissue_mapping_df = fantom5_tissue_mapping_df[['colname','GTEx tissue']]
# Let sample name be the key
fantom5_tissue_mapping_dict = temp_fantom5_tissue_mapping_df.set_index('colname').T.to_dict('list')

i = 0
# Go through each file name in the input directory. Each file name corresponds to a sample
for sample_name in os.listdir(args.output_with_sample_names_path):
    # Take the tissue mapping for each sample
    tissue = fantom5_tissue_mapping_dict[sample_name][0]
    # If tissue mapping isn't na:
    if(not pd.isna(tissue)):
        for t in tissue.split(";"):
            # Save the file with the tissue mapping now appended to the start of the file
            shutil.copyfile(os.path.join(args.output_with_sample_names_path, sample_name), os.path.join(args.output_path, f"{t}_{sample_name}"))
            # Keep track of how many samples were mapped to a tissue correctly
            i = i + 1
            print(i)

print("done")
