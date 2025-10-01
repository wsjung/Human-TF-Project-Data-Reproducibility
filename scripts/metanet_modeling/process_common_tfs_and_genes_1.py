import pandas as pd
import numpy as np
from venny4py.venny4py import *
import os
import argparse

parser = argparse.ArgumentParser(description="processes data and removes autoregulation edges")
parser.add_argument("--tissue", required=True, help="tissue name")
parser.add_argument("--path", required=True, help="input path")
parser.add_argument("--binding_threshold", type=float, default=0.10, help="binding threhsold (default 0.10)")
args = parser.parse_args()

### parameters
tissue = args.tissue
binding_threshold = args.binding_threshold
input_path = args.path

#############
### paths ###
#############
shared_data_path = os.path.join(input_path, "raw_data_shared")
data_path = os.path.join(input_path, "raw_data", tissue)
output_path = os.path.join(input_path, "input_data", tissue, f"{binding_threshold*100}_threshold_metanet")
plot_path = os.path.join(input_path, "plots", tissue)

if not os.path.exists(output_path):
    os.makedirs(output_path)
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

#################
### load data ###
#################
print("loading data")
df_binding = pd.read_csv(os.path.join(shared_data_path, "Remap_binding_scores.bed"), sep="\t")
df_marbach = pd.read_csv(os.path.join(data_path, "MARBACH.tsv"), sep="\t")
# gtex all tissue BART and LASSO
df_bart_all_tissues = pd.read_csv(os.path.join(shared_data_path, "BART_all_gtex_tissues_wide.tsv"), sep="\t", index_col=0)
df_lasso_all_tissues = pd.read_csv(os.path.join(shared_data_path, "LASSO_all_gtex_tissues_wide.tsv"), sep="\t", index_col=0)
df_lasso_all_tissues = df_lasso_all_tissues.T
# gtex tissue-specific data
df_bart_gtex = pd.read_csv(os.path.join(data_path, "net_bart.tsv"), sep="\t", index_col=0)
df_lasso_gtex = pd.read_csv(os.path.join(data_path, "LASSO.txt"), sep="\t", index_col=0)
df_lasso_gtex = df_lasso_gtex.T

###########################
### format/reshape data ### columns: TF, gene, <feature name>
###########################
df_binding.columns = ['TF','gene','binding_score']
df_marbach.columns = ['TF','gene','marbach']
df_bart_all_tissues.index.name='TF'
df_bart_all_tissues = df_bart_all_tissues.reset_index().melt(id_vars='TF', var_name='gene', value_name='bart_all_tissues')
df_bart_gtex.index.name='TF'
df_bart_gtex = df_bart_gtex.reset_index().melt(id_vars='TF', var_name='gene', value_name='bart_gtex')
df_lasso_all_tissues.index.name='TF'
df_lasso_all_tissues = df_lasso_all_tissues.reset_index().melt(id_vars='TF', var_name='gene', value_name='lasso_all_tissues')
df_lasso_gtex.index.name='TF'
df_lasso_gtex = df_lasso_gtex.reset_index().melt(id_vars='TF', var_name='gene', value_name='lasso_gtex')

#######################
### take abs values ###
#######################
print("take abs val")
def take_abs(df):
    df[df.select_dtypes(include=['number']).columns] = df.select_dtypes(include=['number']).abs()
    return df
df_binding = take_abs(df_binding)
df_marbach = take_abs(df_marbach)
df_bart_all_tissues= take_abs(df_bart_all_tissues)
df_bart_gtex = take_abs(df_bart_gtex)
df_lasso_all_tissues= take_abs(df_lasso_all_tissues)
df_lasso_gtex = take_abs(df_lasso_gtex)

###################################
### remove autoregulation edges ###
###################################
print("remove autoregulation edges")
num_autoreg_binding = len(df_binding[(df_binding['TF'] == df_binding['gene']) & (df_binding['binding_score'] > 0)])
num_autoreg_marbach = len(df_marbach[(df_marbach['TF'] == df_marbach['gene']) & (df_marbach['marbach'] > 0)])
num_autoreg_bart_gtex = len(df_bart_gtex[(df_bart_gtex['TF'] == df_bart_gtex['gene']) & (df_bart_gtex['bart_gtex'] > 0)])
num_autoreg_bart_all_tissues = len(df_bart_all_tissues[(df_bart_all_tissues['TF'] == df_bart_all_tissues['gene']) & (df_bart_all_tissues['bart_all_tissues'] > 0)])
num_autoreg_lasso_gtex = len(df_lasso_gtex[(df_lasso_gtex['TF'] == df_lasso_gtex['gene']) & (df_lasso_gtex['lasso_gtex'] > 0)])
num_autoreg_lasso_all_tissues = len(df_lasso_all_tissues[(df_lasso_all_tissues['TF'] == df_lasso_all_tissues['gene']) & (df_lasso_all_tissues['lasso_all_tissues'] > 0)])

df_binding[df_binding['TF'] == df_binding['gene']] = 0
df_marbach[df_marbach['TF'] == df_marbach['gene']] = 0
df_bart_gtex[df_bart_gtex['TF'] == df_bart_gtex['gene']] = 0
df_bart_all_tissues[df_bart_all_tissues['TF'] == df_bart_all_tissues['gene']] = 0
df_lasso_gtex[df_lasso_gtex['TF'] == df_lasso_gtex['gene']] = 0
df_lasso_all_tissues[df_lasso_all_tissues['TF'] == df_lasso_all_tissues['gene']] = 0

#################
### count tfs ###
#################
tfs_binding = np.array(df_binding['TF'].unique(), dtype=str)
tfs_marbach = np.array(df_marbach['TF'].unique(), dtype=str)
tfs_bart_gtex = np.array(df_bart_gtex['TF'].unique(), dtype=str)
tfs_bart_all_tissues = np.array(df_bart_all_tissues['TF'].unique(), dtype=str)
tfs_lasso_gtex = np.array(df_lasso_gtex['TF'].unique(), dtype=str)
tfs_lasso_all_tissues = np.array(df_lasso_all_tissues['TF'].unique(), dtype=str)

###################
### count genes ###
###################
# binding genes are fed in separately since not all genes are guaranteed in the binding BED file
genes_binding = np.loadtxt(os.path.join(shared_data_path, "fantom5_genes_universe_ensg.txt"), dtype=str)
genes_marbach = np.array(df_marbach['gene'].unique(), dtype=str)
genes_bart_gtex = np.array(df_bart_gtex['gene'].unique(), dtype=str)
genes_bart_all_tissues = np.array(df_bart_all_tissues['gene'].unique(), dtype=str)
genes_lasso_gtex = np.array(df_lasso_gtex['gene'].unique(), dtype=str)
genes_lasso_all_tissues = np.array(df_lasso_all_tissues['gene'].unique(), dtype=str)
print("## TF/gene universe by datset")
print(f"binding: {len(tfs_binding)} TFs x {len(genes_binding)} gene")
print(f"marbach: {len(tfs_marbach)} TFs x {len(genes_marbach)} gene")
print(f"bart_gtex: {len(tfs_bart_gtex)} TFs x {len(genes_bart_gtex)} gene")
print(f"bart_all_tissues: {len(tfs_bart_all_tissues)} TFs x {len(genes_bart_all_tissues)} gene")
print(f"lasso_gtex: {len(tfs_lasso_gtex)} TFs x {len(genes_lasso_gtex)} gene")
print(f"lasso_all_tissues: {len(tfs_lasso_all_tissues)} TFs x {len(genes_lasso_all_tissues)} gene")

#################################
### find common tfs and genes ###
#################################
common_tfs = set(tfs_binding).intersection(tfs_marbach).intersection(tfs_bart_gtex)
common_genes = set(genes_binding).intersection(genes_marbach).intersection(genes_bart_gtex)
print(f"common: {len(common_tfs)} TFs x {len(common_genes)} genes")

#########################
### plot venn diagram ###
#########################
print("plotting venn diagrams")
sets_tfs = {
    'binding': set(tfs_binding),
    'marbach': set(tfs_marbach),
    'GTEx (BART/LASSO)': set(tfs_bart_gtex)
}
sets_genes = {
    'binding': set(genes_binding),
    'marbach': set(genes_marbach),
    'GTEx (BART/LASSO)': set(genes_bart_gtex)
}
venn_path_tf = os.path.join(plot_path, "venn_diagrams", "tfs")
venn_path_gene = os.path.join(plot_path, "venn_diagrams", "genes")
if not os.path.exists(venn_path_tf):
    os.makedirs(venn_path_tf)
if not os.path.exists(venn_path_gene):
    os.makedirs(venn_path_gene)
venny4py(sets=sets_tfs, out=venn_path_tf)
venny4py(sets=sets_genes, out=venn_path_gene)

########################################
### subset data for common tfs/genes ###
########################################
print("subsetting for common tfs/genes")
df_binding_subset = df_binding[(df_binding['TF'].isin(common_tfs)) & (df_binding['gene'].isin(common_genes))]
df_marbach_subset = df_marbach[(df_marbach['TF'].isin(common_tfs)) & (df_marbach['gene'].isin(common_genes))]
df_bart_gtex_subset = df_bart_gtex[(df_bart_gtex['TF'].isin(common_tfs)) & (df_bart_gtex['gene'].isin(common_genes))]
df_bart_all_tissues_subset = df_bart_all_tissues[(df_bart_all_tissues['TF'].isin(common_tfs)) & (df_bart_all_tissues['gene'].isin(common_genes))]
df_lasso_gtex_subset = df_lasso_gtex[(df_lasso_gtex['TF'].isin(common_tfs)) & (df_lasso_gtex['gene'].isin(common_genes))]
df_lasso_all_tissues_subset = df_lasso_all_tissues[(df_lasso_all_tissues['TF'].isin(common_tfs)) & (df_lasso_all_tissues['gene'].isin(common_genes))]

#############
### merge ###
#############
print("merging")
df_merged_subset = df_binding_subset.merge(
                df_marbach_subset, on=['TF','gene'], how='outer').merge(
                df_bart_gtex_subset, on=['TF','gene'], how='outer').merge(
                df_bart_all_tissues_subset, on=['TF','gene'], how='outer').merge(
                df_lasso_gtex_subset, on=['TF','gene'], how='outer').merge(
                df_lasso_all_tissues_subset, on=['TF','gene'], how='outer')
# fill any NA's with 0's
df_merged_subset = df_merged_subset.fillna(0)

#############################
### filter binding labels ###
#############################
print(f"filtering binding labels for max {binding_threshold*100}% positive labels")
num_tfs = len(common_tfs)
num_genes = len(common_genes)
num_genes_threshold = int(num_genes * binding_threshold)

df_label_list = []
num_pos_prefilter = []
num_pos_postfilter = []
threshold_count = 0
for tf in df_merged_subset['TF'].unique():
    df_tf = df_merged_subset[df_merged_subset['TF'] == tf]
    df_tf['label'] = np.where(df_tf['binding_score'] > 0, 1, 0)

    # num positive labels before filtering
    num_pos = df_tf['label'].sum()
    num_pos_prefilter.append(num_pos)

    if num_pos > num_genes_threshold:
        threshold_count = threshold_count + 1
        df_tf = df_tf.sort_values(by='binding_score', ascending=False)
        df_tf.iloc[num_genes_threshold:, df_tf.columns.get_loc('label')] = 0

    # num positive labels after filtering
    num_pos = df_tf['label'].sum()
    num_pos_postfilter.append(num_pos)

    df_label_list.append(df_tf)
df_label = pd.concat(df_label_list)
print(f"{threshold_count} TFs had more positive labels than {binding_threshold*100}% of gene universe")

plt.figure(figsize=(10,5))
plt.hist(num_pos_prefilter, bins=30, label='prefilter')
plt.hist(num_pos_postfilter, bins=30, label='postfilter')
plt.legend()
plt.title(f"{tissue}\nnumber of positive labels pre/post per-TF {binding_threshold*100}% filtering ({num_tfs} TFs x {num_genes} genes)")
plt.xlabel("num positive labels")
plt.ylabel("num TFs")
plt.savefig(os.path.join(plot_path, f"num_labels_per_tf_{binding_threshold*100}_threshold_metanet.png"))

##################################
### rename cols with uppercase ###
##################################
df_label.columns = map(str.upper, df_label.columns)

#####################
### write to file ###
#####################
df_label.to_csv(os.path.join(output_path, f"{tissue}_{binding_threshold*100}_threshold_metanet.txt"), sep="\t", index=False)


print("DONE")
