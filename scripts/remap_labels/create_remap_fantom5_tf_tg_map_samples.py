import pandas as pd
import numpy as np
import os
import re
import argparse

parser = argparse.ArgumentParser(description="maps REMAP TF->promoter and promoter-gene networks to TF->gene network")
parser.add_argument("--gencode_protein", required=True, help="file with gencode protein coding genes")
parser.add_argument("--REMAP_promoter", required=True, help="REMAP TF peaks mapped to promoter regions")
parser.add_argument("--promoter_gene", required=True, help="file with FANTOM5 promoters mapped to genes")
parser.add_argument("--gtex_tissue_mapping", required=True, help="file with FANTOM5 samples to GTEx tissue mapping")
parser.add_argument("--output_dir", required=True, help="output directory")
args = parser.parse_args()

OUTPUT_DIR = os.path.join(args.output_dir, "tissues")
OUTPUT_DIR_GENCODE_ID_ZERO_FILTER = os.path.join(args.output_dir, "tissues_tf_gencode") # TFs with gene_id and no filtering
OUTPUT_DIR_GENCODE_ID_75_FILTER = os.path.join(args.output_dir, "tissues_tf_gencode_75percentile") # TFs with gene_id and edges filtered for top 25 percentile

# load gene name mapping
df_gencode29_protein_coding_genes = pd.read_csv(args.gencode_protein, sep="\t")
df_gencode29_dict = {row['gene_name']:row['gene_id'] for _,row in df_gencode29_protein_coding_genes.iterrows()}

# load TF to promoter data
df_remap_tf_to_promoter = pd.read_csv(args.REMAP_promoter, sep="\t")

# load promoter to gene data
## promoters mapped to gene_ids
## first two rows (stats rows) dropped
df_promoter_to_gene = pd.read_csv(args.promoter_gene, sep="\t")

# load GTEx tissue mapping
df_tissue_mapping = pd.read_csv(args.gtex_tissue_mapping, sep="\t")
## explode GTEx tissue: some samples are mapped to multiple GTEx tissues
df_tissue_mapping['GTEx tissue'] = df_tissue_mapping['GTEx tissue'].str.split(';')
df_tissue_mapping = df_tissue_mapping.explode('GTEx tissue')

fantom5_samples = df_tissue_mapping['colname'].unique()
gtex_tissues = df_tissue_mapping['GTEx tissue'].dropna().unique()
print(f"{len(fantom5_samples)} unique FANTOM5 samples")
print(f"{len(gtex_tissues)} unique GTEx tissues")

# for each gtex tissue
for gtex_tissue in gtex_tissues:
    print(f"##########\n### GTEx tissue: {gtex_tissue}\n##########")

    # replace special chars with underscores
    gtex_tissue_underscore = re.sub(r'[^\w.]', '_', gtex_tissue)

    # FANTOM5 samples that map to gtex tissue
    df_gtex_tissue_samples = df_tissue_mapping[df_tissue_mapping['GTEx tissue'] == gtex_tissue]
    gtex_tissue_samples = df_gtex_tissue_samples['colname'].unique()
    print(f"{len(gtex_tissue_samples)} samples map to {gtex_tissue}")

    # subset promoter->gene by samples
    subset_cols = ['00Annotation','gene_id'] + list(gtex_tissue_samples)
    subset_cols = set(subset_cols).intersection(df_promoter_to_gene.columns) # subset_cols colnames must exist in df
    df_promoter_to_gene_tissue = df_promoter_to_gene[subset_cols]
    #print(df_promoter_to_gene)
    #print(df_promoter_to_gene_tissue)

    # merge with TF->promoter data
    df_tf_to_gene = df_remap_tf_to_promoter.merge(df_promoter_to_gene_tissue, left_on='name', right_on='00Annotation', how='inner')
    df_tf_to_gene = df_tf_to_gene.drop(columns=['chr','start','end','name','00Annotation'])

    print("filtering")
    # sum activity levels across samples
    df_tf_to_gene['sum'] = df_tf_to_gene[gtex_tissue_samples].sum(axis=1)

    # TF to gencode mapping
    print("mapping TFs to gencode gene_id")
    df_tf_to_gene['TF_gene_id'] = df_tf_to_gene['TF'].str.upper().map(df_gencode29_dict)
    num_TFs_lost = df_tf_to_gene['TF'].nunique() - df_tf_to_gene['TF_gene_id'].nunique()
    print(f"> lost {num_TFs_lost} TFs")

    ## only keep TF-TG pairs where sum of normalized activity > 0
    df_tf_to_gene_zero_filtering = df_tf_to_gene[df_tf_to_gene['sum'] > 0]

    ## only keep TF-TG pairs where sum of norm activity >= 75 percentile
    _75_percentile_filtering_threshold = df_tf_to_gene['sum'].quantile(.75)
    df_tf_to_gene_75_filtering = df_tf_to_gene[df_tf_to_gene['sum'] > _75_percentile_filtering_threshold]

    # only write TF and TG (zero filtering)
    df_tf_to_gene_zero_filtering = df_tf_to_gene_zero_filtering[['TF_gene_id','gene_id']].drop_duplicates()
    df_tf_to_gene_zero_filtering = df_tf_to_gene_zero_filtering.rename(columns={'TF_gene_id':'regulator','gene_id':'target'})

    # only write TF and TG (75% filtering)
    df_tf_to_gene_75_filtering = df_tf_to_gene_75_filtering[['TF_gene_id','gene_id']].drop_duplicates()
    df_tf_to_gene_75_filtering = df_tf_to_gene_75_filtering.rename(columns={'TF_gene_id':'regulator','gene_id':'target'})

    print("mapping gene names to gene_ids")
    # map TF to gene_id (zero filtering)
    df_tf_to_gene_zero_filtering_ensg = df_tf_to_gene_zero_filtering
    #df_tf_to_gene_zero_filtering_ensg['regulator'] = df_tf_to_gene_zero_filtering['regulator'].str.upper()map(df_gencode29_dict)
    df_tf_to_gene_zero_filtering_ensg = df_tf_to_gene_zero_filtering_ensg.dropna()

    # map TF to gene_id (75% filtering)
    df_tf_to_gene_75_filtering_ensg = df_tf_to_gene_75_filtering
    #df_tf_to_gene_75_filtering_ensg['regulator'] = df_tf_to_gene_75_filtering['regulator'].str.upper().map(df_gencode29_dict)
    df_tf_to_gene_75_filtering_ensg = df_tf_to_gene_75_filtering_ensg.dropna()

    print("writing to disk")
    # write to disk
    df_tf_to_gene_zero_filtering.to_csv(os.path.join(OUTPUT_DIR, f"{gtex_tissue_underscore}.txt"), sep="\t", index=False)

    # write to disk (tfs with gencode29 gene_id)
    df_tf_to_gene_zero_filtering_ensg.to_csv(os.path.join(OUTPUT_DIR_GENCODE_ID_ZERO_FILTER, f"{gtex_tissue_underscore}.txt"), sep="\t", index=False)

    # write to disk (75% percentile filtering) + (tfs with gencode29 gene_id)
    df_tf_to_gene_75_filtering_ensg.to_csv(os.path.join(OUTPUT_DIR_GENCODE_ID_75_FILTER, f"{gtex_tissue_underscore}.txt"), sep="\t", index=False)






print("DONE")
