import pandas as pd
import numpy as np
import itertools
import os
import re
import argparse

parser = argparse.ArgumentParser(description="creates TF->TG matrix based on enhanceratlas FIMO scores")
parser.add_argument("--gtex_mapping", required=True, help="GTEx tissue mapping file for fantom5 samples")
parser.add_argument("--ea_name_mapping", required=True, help="mapping file for enhanceratlas hg19 and hg38 names")
parser.add_argument("--gencode_id_mapping", required=True, help="mapping file between gencode protein-coding gene name, IDs")
parser.add_argument("--fimo", required=True, help="normalized fimo scores")
parser.add_argument("--ea_gene_interactions", required=True, help="normalized enhanceratlas gene interaction scores")
parser.add_argument("--output_dir", required=True, help="output directory")
parser.add_argument("--output_dir_gencode", required=True, help="output directory for matrix with gencode TF gene_ids")
args = parser.parse_args()


# mapping files
df_GTEx_mapping = pd.read_csv(args.gtex_mapping, sep="\t")
df_sequence_name_mapping = pd.read_csv(args.ea_name_mapping, header=None, names=['sequence_name_hg19','sequence_name_hg38'])
df_gencode29_protein_coding_genes = pd.read_csv(args.gencode_id_mapping, sep="\t")
df_gencode29_dict = {row['gene_name']:row['gene_id'] for _,row in df_gencode29_protein_coding_genes.iterrows()}

gtex_tissues = df_GTEx_mapping['GTEx'].unique()
print(f"# {len(gtex_tissues)} unique GTEx tissues")

# data
df_fimo_all = pd.read_csv(args.fimo, sep="\t") # normalized
df_ge_all = pd.read_csv(args.ea_gene_interactions, sep="\t") # normalized
## HeLa-S3 is incorrectly labelled as Hela-S3
df_ge_all['cell'] = df_ge_all['cell'].replace({'Hela-S3':'HeLa-S3'})

## map hg19 sequencename for fimo results
df_fimo_all = df_fimo_all.merge(df_sequence_name_mapping, on='sequence_name_hg38', how='inner')


# map gene interaction genes to gencode ENSG ids
## use gene_id
df_ge_all = df_ge_all.merge(df_gencode29_protein_coding_genes, left_on='gene_ensg', right_on='ENSG_ID_enhanceratlas', how='inner')

# create an empty df of all TFxTG pairs
tfs_all = set(df_fimo_all['motif_alt_id'])
print(f"### {len(tfs_all)} total TFs")
genes_all = set(df_gencode29_protein_coding_genes['gene_id'])
df_zero_edges = pd.DataFrame(list(itertools.product(tfs_all,genes_all)), columns=['motif_alt_id','gene_id'])
df_zero_edges['tf_tg_score'] = 0


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

    # subset fimo hits in cells that map to gtex tissue
    df_fimo_tissue = df_fimo_all[df_fimo_all['cell'].isin(gtex_tissue_cells)]
    fimo_missing_cells = set(gtex_tissue_cells) - set(df_fimo_all['cell'].unique())
    if len(fimo_missing_cells) > 0:
        print(f"# {len(fimo_missing_cells)} cell(s) missing from FIMO results:")
        print(list(fimo_missing_cells))

    # subset gene interactions
    df_ge_tissue = df_ge_all[df_ge_all['cell'].isin(gtex_tissue_cells)]
    ge_missing_cells = set(gtex_tissue_cells) - set(df_ge_tissue['cell'].unique())
    if (len(ge_missing_cells) > 0):
        print(f"# {len(ge_missing_cells)} cell(s) missing from gene interaction scores:")
        print(list(ge_missing_cells))

    # merge
    df_tissue = df_fimo_tissue.merge(df_ge_tissue, on='sequence_name_hg19', how='inner', suffixes=('_fimo','_interaction'))

    # calculate tf->tg score (motif_score_norm x interaction_score_norm)
    df_tissue['tf_tg_score'] = df_tissue['motif_score_norm'] * df_tissue['interaction_score_norm']

    #print(f"num TFs: {df_tissue['motif_alt_id'].nunique()}")
    #print(f"num TGs: {df_tissue['gene_id'].nunique()}")

    df_tissue = df_tissue[['motif_alt_id','gene_id','motif_score_norm','interaction_score_norm','tf_tg_score','sequence_name_hg19','cell_fimo']]
    #print(df_tissue)

    #print(df_tissue[df_tissue[['motif_alt_id','gene_id']].duplicated(keep=False)].sort_values(by=['motif_alt_id','gene_id']))

    # take maximum score for each pair of TF-TG
    df_tissue = df_tissue.groupby(['motif_alt_id','gene_id'], as_index=False)['tf_tg_score'].max()
    #print(df_tissue)

    # fill in zero edges
    df_tissue = pd.concat([df_tissue, df_zero_edges])
    df_tissue = df_tissue.groupby(['motif_alt_id','gene_id'], as_index=False)['tf_tg_score'].max()
    #print(df_tissue)

    # convert tf name to gene_id
    df_tissue['motif_alt_id'] = df_tissue['motif_alt_id'].str.upper()
    df_tissue['motif_alt_id_gene_id'] = df_tissue['motif_alt_id'].map(df_gencode29_dict)
    ## number of TFs lost to mapping:
    #print(f"tfs lost to mapping: {df_tissue[df_tissue['motif_alt_id_gene_id'].isna()]['motif_alt_id'].unique()}")

    # add in TFs with 0 scores


    # long to wide
    df_tissue_wide = df_tissue.pivot(columns='gene_id', index='motif_alt_id', values='tf_tg_score')
    # long to wide (tfs with gencode29 gene_id)
    df_tissue_wide_tf_gencode_id = df_tissue.pivot(columns='gene_id', index='motif_alt_id_gene_id', values='tf_tg_score')

    # write to disk
    df_tissue_wide.to_csv(os.path.join(args.output_dir, f"{gtex_tissue_underscore}.txt"), sep="\t", index=True)
    # write to disk (tfs with gencode29 gene_id)
    df_tissue_wide_tf_gencode_id.to_csv(os.path.join(args.output_dir_gencode, f"{gtex_tissue_underscore}.txt"), sep="\t", index=True)


print(f"DONE")


