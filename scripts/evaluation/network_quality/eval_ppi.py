import pandas as pd
import numpy as np
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(description="calculates PPI support")
parser.add_argument("--input_net", required=True, help="path to input network, 3 tab-separated cols for TFs, genes, scores")
parser.add_argument("--PPI_db", required=True, help="path to PPI database to use as evaluation gold-standard, space-separated")
parser.add_argument("--top_TF_pairs", default=100, type=int, help="number of top (TF,TF) pairs to evaluate (default=100)")
parser.add_argument("--pairs_per_threshold", default=10, type=int, help="number of (TF,TF) pairs to evaluate per threshold (default=10)")
parser.add_argument("--genes_per_TF", default=250, type=int, help="number of top edges to evaluate per TF")
parser.add_argument("--output_path", required=True, help="output path")
args = parser.parse_args()


def calc_jaccard_sim(set1, set2):

    intersection = set(set1).intersection(set(set2))
    union = set(set1).union(set(set2))
    jaccard_sim = len(intersection) / len(union) if union else 0
    return jaccard_sim


#################
# read networks #
#################
df_network = pd.read_csv(args.input_net, sep="\t", header=None, names=['TF','GENE','SCORE'])
num_tfs = df_network['TF'].nunique()
num_genes = df_network['GENE'].nunique()

df_ppi = pd.read_csv(args.PPI_db, sep=" ")

##############
# preprocess #
##############
## PPI
# ensure stable ID
df_ppi['protein1'] = df_ppi['protein1'].str.split('.', expand=True)[0]
df_ppi['protein2'] = df_ppi['protein2'].str.split('.', expand=True)[0]

# subset confidence score >= 0.7
df_ppi = df_ppi[df_ppi['combined_score'] >= 700]

# set of all TFs in PPI db
ppi_tfs = set(df_ppi['protein1'].unique()).union(df_ppi['protein2'].unique())
print(f" - {len(ppi_tfs)} proteins found in PPI db")

## network
# subset for TFs in PPI database
num_tfs = df_network['TF'].nunique()
num_genes = df_network['GENE'].nunique()
print(f" - {num_tfs} TFs found in input network")
print(" > subsetting for TFs in PPI db")
df_network = df_network[df_network['TF'].isin(ppi_tfs)]
num_tfs = df_network['TF'].nunique()
print(f" - {num_tfs} TFs in input network found in PPI db")

# subset top num_tfs * genes_per_TF edges
df_network['SCORE'] = df_network['SCORE'].abs()
df_network = df_network.sort_values(by='SCORE', ascending=False)
num_top_edges = num_tfs * args.genes_per_TF
top_edges_idx = min(num_top_edges, len(df_network)) # whichever's smaller of total number of edges and top edges threshold
df_network = df_network.iloc[:top_edges_idx,]

################################
# calculate jaccard similarity #
################################
print("calculating jaccard similarity")
# build dictionary mapping TF to set of GENEs
tf_to_genes_dict = df_network.groupby('TF')['GENE'].apply(set).to_dict()

# get unique (TF,TF) pairs
tfs = list(tf_to_genes_dict.keys())
tf_pairs = combinations(tfs,2)

results = []
for tf1, tf2 in tf_pairs:
    genes_tf1 = tf_to_genes_dict[tf1]
    genes_tf2 = tf_to_genes_dict[tf2]
    jaccard_sim = calc_jaccard_sim(genes_tf1, genes_tf2)
    results.append({
        'TF1': tf1,
        'TF2': tf2,
        'jaccard': jaccard_sim
    })
df_jaccard = pd.DataFrame(results)
# sort by jaccard similarity
df_jaccard = df_jaccard.sort_values(by='jaccard', ascending=False)
print(df_jaccard)


#####################
# calculate support #
#####################
print("calculating PPI support")
# set of PPI support pairs
ppi_pairs = set(
    frozenset([ row['protein1'], row['protein2'] ])
    for _, row in df_ppi.iterrows()
)

# function to check if TF pair is in ppi_pairs
def _is_supported_by_ppi(row):
    tf_pair = frozenset([ row['TF1'], row['TF2'] ])
    return tf_pair in ppi_pairs

# apply check
df_jaccard['supported_by_ppi'] = df_jaccard.apply(_is_supported_by_ppi, axis=1)
print(df_jaccard)

##################################
# calculate support by threshold #
##################################
# sort by jaccard similarity
df_jaccard = df_jaccard.sort_values(by='jaccard', ascending=False).reset_index(drop=True)

# convert support to int
df_jaccard['supported_by_ppi_int'] = df_jaccard['supported_by_ppi'].astype(int)

# calculate cumulative sum
df_jaccard['cumulative_ppi_support'] = df_jaccard['supported_by_ppi_int'].cumsum()
df_jaccard['cumulative_count'] = df_jaccard.index + 1

# caulcate cumulative PPI support ratio
df_jaccard['cumulative_ppi_support_percentage'] = df_jaccard['cumulative_ppi_support'] / df_jaccard['cumulative_count'] * 100

# subset thresholds
threshold_idx = [i for i in range(9, args.top_TF_pairs, args.pairs_per_threshold)]

# extract cumulative PPI support percentages
df_jaccard_final = df_jaccard.loc[threshold_idx, ['cumulative_count','cumulative_ppi_support_percentage']].reset_index(drop=True)
df_jaccard_final.columns = ['rank','ppi_support']
print(df_jaccard_final)

# write to disk
df_jaccard_final.to_csv(args.output_path, sep="\t", index=False)

print("DONE")
