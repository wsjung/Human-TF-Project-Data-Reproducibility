import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Convert genie3 output into correct matrices')
parser.add_argument('--genie3_output', required=True, help='Path to the genie3 outputs for each gene.')
parser.add_argument('--genes_list', required=True, help='newline separated list of all the gene ensg ids')
parser.add_argument('--tfs_list', required=True, help='all the tfs in the project')
parser.add_argument('--output_dir', required=True, help='Directory to save the output CSV files.')
args = parser.parse_args()

genes_list_df = pd.read_csv(args.genes_list, sep='\t', header=None)
tfs_list_df = pd.read_csv(args.tfs_list, sep='\t')

genes_list = genes_list_df[0]

zeros_matrix = np.zeros((len(genes_list),len(genes_list)))
matrix_df = pd.DataFrame(zeros_matrix, columns = genes_list, index = genes_list)

def add_to_matrix(gene):
    gene_output_df = pd.read_csv(os.path.join(args.genie3_output, f'{gene}.csv'), sep=',')
    matrix_df.loc[gene, :] = gene_output_df['feature_importance'].values

genes_list.apply(add_to_matrix)

matrix_df.index.name = None

matrix_df.to_csv(os.path.join(args.output_dir, f"gene_by_gene.csv"), sep=",", index=True)

tfs_matrix_df = matrix_df[matrix_df.index.isin(tfs_list_df['TF_ensg'])]
tfs_matrix_df.to_csv(os.path.join(args.output_dir, f"tf_by_gene.csv"), sep=",", index=True)

print('done')
