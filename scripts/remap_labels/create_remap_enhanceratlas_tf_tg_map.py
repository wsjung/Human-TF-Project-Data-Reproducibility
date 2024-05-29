import pandas as pd
import numpy as np
import itertools
import os
import re
import argparse

parser = argparse.ArgumentParser(description="creates REMAP TF->TG labels via enhancer atlas enhancer regions")
parser.add_argument("--enhanceratlas_gtex_mapping", required=True, help="mapping file for enhanceratlas tissues to GTEx tissues")
parser.add_argument("--gene_interaction_gencode_mapping", required=True, help="gencode mapping for enhanceratlas gene interactions")
parser.add_argument("--sequence_name_mapping", required=True, help="mapping between hg19 and lifted hg38 names for enhanceratlas regions")
parser.add_argument("--gencode_protein", required=True, help="gencode protein coding genes")
parser.add_argument("--gene_interactions", required=True, help="enhanceratlas gene interactions (normalized)")
parser.add_argument("--REMAP_directory", required=True, help="directory containing REMAP TF peaks to enhanceratlas regions named by celltype")
parser.add_argument("--output_dir", required=True, help="output directory")
args = parser.parse_args()


OUTPUT_DIR = os.path.join(args.output_dir, "tissues")
OUTPUT_DIR_GENCODE_ID = os.path.join(args.output_dir, "tissues_tf_gencode") # TFs with gencode gene_id

# mapping files
df_GTEx_mapping = pd.read_csv(args.enhanceratlas_gtex_mapping, sep="\t")
df_ge_gencode_mapping = pd.read_csv(args.gene_interaction_gencode_mapping, sep="\t")
gene_name_id_dict = {row['gene_name']:row['gene_id'] for _,row in df_ge_gencode_mapping.iterrows()}
df_sequence_name_mapping = pd.read_csv(args.sequence_name_mapping, header=None, names=['sequence_name_hg19','sequence_name_hg38'])
df_gencode29_protein_coding_genes = pd.read_csv(args.gencode_protein, sep="\t")
df_gencode29_dict = {row['gene_name']:row['gene_id'] for _,row in df_gencode29_protein_coding_genes.iterrows()}

gtex_tissues = df_GTEx_mapping['GTEx'].unique()
print(f"# {len(gtex_tissues)} unique GTEx tissues")


# data
### gene interactions
df_ge_all = pd.read_csv(args.gene_interactions, sep="\t") # normalized
## HeLa-S3 is incorrectly labelled as Hela-S3
df_ge_all['cell'] = df_ge_all['cell'].replace({'Hela-S3':'HeLa-S3'})
df_ge_all['sequence_name_hg19'] = df_ge_all['sequence_name_hg19'].str.replace("Hela_S3","HeLa-S3")
### remap
remap_by_cell_dir = args.REMAP_directory
remap_cells = os.listdir(remap_by_cell_dir)

# map gene interaction genes to gencode ENSG ids
## use gene_id
df_ge_all = df_ge_all.merge(df_ge_gencode_mapping, left_on='gene_ensg', right_on='ENSG_ID_enhanceratlas', how='inner')

# for each gtex tissue
for gtex_tissue in gtex_tissues:

    print(f"### GTEx tissue: {gtex_tissue} ###")

    # replace tissue names' special characters with underscores
    gtex_tissue_underscore = re.sub(r'[^\w.]', '_', gtex_tissue)

    # subset cells that map to gtex tissue
    df_tissue_mapping = df_GTEx_mapping[df_GTEx_mapping['GTEx'] == gtex_tissue]
    gtex_tissue_cells = df_tissue_mapping['Cells'].unique()
    #print(f"# {len(gtex_tissue_cells)} enhancer atlas cell(s) map to {gtex_tissue}")
    print(list(gtex_tissue_cells))

    # subset gene interactions
    df_ge_tissue = df_ge_all[df_ge_all['cell'].isin(gtex_tissue_cells)]
    ge_missing_cells = set(gtex_tissue_cells) - set(df_ge_tissue['cell'].unique())
    if (len(ge_missing_cells) > 0):
        print(f"# {len(ge_missing_cells)} cell(s) missing from gene interaction scores:")
        print(list(ge_missing_cells))


    # find remap hits in gtex tissue cells
    remap_missing_cells = set(gtex_tissue_cells) - set(remap_cells)
    remap_common_cells = set(gtex_tissue_cells).intersection(remap_cells)
    df_remap = []
    for remap_cell in remap_common_cells:
        df = pd.read_csv(os.path.join(remap_by_cell_dir, remap_cell, f"remap_enhancers.{remap_cell}.bed"), sep="\t")
        df['cell'] = remap_cell
        df_remap.append(df)
    df_remap = pd.concat(df_remap)
    df_remap = df_remap.sort_values(by=['chr','start','end'])

    if len(remap_missing_cells) > 0:
        print(f"# {len(remap_missing_cells)} cell(s) missing from FIMO results:")
        print(list(remap_missing_cells))

    # merge
    df_tissue = df_remap.merge(df_ge_tissue, left_on='name', right_on='sequence_name_hg19', how='inner', suffixes=('_remap','_interaction'))
    df_tissue = df_tissue[['TF','gene_id','interaction_score_norm','sequence_name_hg19','cell_remap']]

    # take maximum score for each pair of TF-TG
    df_tissue = df_tissue.groupby(['TF','gene_id'], as_index=False)['interaction_score_norm'].max()
    df_tissue['TF'] = df_tissue['TF'].str.upper()
    df_tissue['TF_gene_id'] = df_tissue['TF'].map(df_gencode29_dict)

    ## number of TFs lost from mapping to gencode IDs
    lost_TFs = df_tissue[df_tissue['TF_gene_id'].isna()]['TF'].unique()
    num_TFs_lost = df_tissue['TF'].nunique() - df_tissue['TF_gene_id'].nunique()
    num_edges_lost = len(df_tissue) - len(df_tissue.dropna(subset=['TF_gene_id']))
    print(f"> lost {num_TFs_lost} TFs ({num_edges_lost})")
    print(f">> {lost_TFs}")

    # write to disk
    df_tissue_tf_hgnc = df_tissue[['TF','gene_id']].rename(columns={'TF':'regulator','gene_id':'target'})
    df_tissue_tf_hgnc.to_csv(os.path.join(OUTPUT_DIR, f"{gtex_tissue_underscore}.txt"), sep="\t", index=False)

    # write to disk (tfs with gencode29 gene_id)
    df_tissue_tf_gencode_id = df_tissue[['TF_gene_id','gene_id']].rename(columns={'TF_gene_id':'regulator', 'gene_id':'target'})
    df_tissue_tf_gencode_id = df_tissue_tf_gencode_id.dropna()
    df_tissue_tf_gencode_id.to_csv(os.path.join(OUTPUT_DIR_GENCODE_ID, f"{gtex_tissue_underscore}.txt"), sep="\t", index=False)










print("DONE")
