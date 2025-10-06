import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="processing of PANDA networks")
parser.add_argument("--panda_path", required=True, help="path containing raw PANDA matrices extracted from .RData file")
parser.add_argument("--metanet_path", required=True, help="path to METANet networks")
parser.add_argument("--gene_mapping", required=True, help="path to gene symbol and ID mapping TSV file")
args = parser.parse_args()

#############
# load data #
#############
## PANDA matrices
df_edges = pd.read_csv(os.path.join(panda_raw_path, "Sonawane_edges.tsv"), sep="\t")
df_gene_annt = pd.read_csv(os.path.join(panda_raw_path, "Sonawane_genes_annotation.tsv"), sep="\t")
df_tissue_annotation = pd.read_csv(os.path.join(panda_raw_path, "Sonawane_network_tissue_specificity_annotation_matrix.tsv"), sep="\t")
df_weights = pd.read_csv(os.path.join(panda_raw_path, "Sonawane_network_weights_matrix.tsv"), sep="\t")

## gene mapping file
df_gencode = pd.read_csv(args.gene_mapping, sep="\t")
df_gencode = df_gencode[~df_gencode['gene_id'].str.contains('PAR_Y')] # remove PAR_Y genes
df_gencode['gene_id_stable'] = df_gencode['gene_id'].str.split('.', expand=True)[0] # stable ID

##############
# Preprocess #
##############
# Map PANDA tissues to METANet / GTEx tissue names
tissues_panda = list(df_tissue_annotation.columns)
panda_xgboost_tissue_mapping_dict = {
    'adipose_visceral':'adipose_visceral_omentum',
    'breast':'breast_mammary_tissue',
    'skeletal_muscle':'muscle_skeletal',
    'tibial_nerve':'nerve_tibial'
}
## lowercase PANDA tissues
tissues_panda_lower = [tissue.lower() for tissue in tissues_panda]
tissues_panda_lower_mapped = [
    panda_xgboost_tissue_mapping_dict[tissue] if tissue in panda_xgboost_tissue_mapping_dict.keys() else tissue 
    for tissue in tissues_panda_lower
]

# subset for protein-coding genes
df_edges_protein_coding = df_edges.copy()
df_edges_protein_coding = df_edges_protein_coding[df_edges_protein_coding['Gene'].isin(df_gencode['gene_id_stable'])]

# map TF HGNC name to TF ENSG stable ID
df_panda_gene_annt_dict = { row['Symbol']:row['Name'] for _,row in df_gene_annt.iterrows() }
df_edges_protein_coding['TF_ENSG'] = df_edges_protein_coding['TF'].map(df_panda_gene_annt_dict)

# remove self-regulation edges
df_edges_protein_coding_noselfreg = df_edges_protein_coding.copy()
df_edges_protein_coding_noselfreg = df_edges_protein_coding_noselfreg[df_edges_protein_coding_noselfreg['TF_ENSG'] != df_edges_protein_coding_noselfreg['Gene']]

############################
# Subset and save networks #
############################

n_tfs_common_list = []
for tissue_panda,tissue_metanet in zip(tissues_panda,tissues_panda_lower_mapped):
    print(f"# {tissue_panda} ({tissue_metanet})")

    # create tissue df from filtered edges
    df_tissue = df_edges_protein_coding_noselfreg.copy()
    
    # extract edge weights for tissue
    tissue_idx = df_tissue.index
    tissue_edge_weights = df_weights.loc[tissue_idx, tissue_panda]

    # append edge weights
    df_tissue['SCORE'] = tissue_edge_weights

    # sort edge weight descending
    df_tissue = df_tissue.sort_values(by='SCORE', ascending=False)

    # subset columns
    df_tissue = df_tissue[['TF_ENSG','Gene','SCORE']]
    df_tissue.columns = ['TF','GENE','SCORE']

    # save
    print("> saving panda network")
    panda_tissue_output_path = os.path.join(args.metanet_path, tissue_metanet, "benchmark_networks", "panda")
    os.makedirs(panda_tissue_output_path, exist_ok=True)
    df_tissue.to_csv(os.path.join(panda_tissue_output_path, f"{tissue_metanet}_panda_stableID.tsv"), sep="\t", index=False, header=False)

    # calculate max targets per TF
    max_targets_per_tf = int(len(df_tissue) / df_tissue['TF'].nunique())
    print(f"\t> max {max_targets_per_tf} targets per TF")

    # TISSUE SPECIFIC NETWORK
    df_tissue_annotation_idx = df_tissue_annotation[df_tissue_annotation[tissue_panda] == 1].index
    df_tissue_specific = df_edges_protein_coding_noselfreg.copy()
    final_tissue_specific_idx = np.intersect1d(tissue_idx, df_tissue_annotation_idx) # intersect tissue specific edge idx and filtered edge idx
    df_tissue_specific = df_tissue_specific.loc[final_tissue_specific_idx]
    df_tissue_specific['SCORE'] = df_weights.loc[final_tissue_specific_idx, tissue_panda]
    df_tissue_specific = df_tissue_specific.sort_values(by='SCORE', ascending=False)
    df_tissue_specific = df_tissue_specific[['TF_ENSG','Gene','SCORE']]
    df_tissue_specific.columns = ['TF','GENE','SCORE']


    # subset for TFs in matching METANet
    df_tissue_subset = df_tissue.copy()
    try:
        df_metanet = pd.read_csv(os.path.join(xgboost_path, tissue_metanet, "METANet", "xgboost", f"{tissue_metanet}_xgboost_stableID.tsv"), sep="\t", header=None, names=['TF','GENE','SCORE'])
        df_metanet_subset = df_metanet.copy()

        TFs_common = set(df_metanet['TF'].unique()).intersection(set(df_tissue_subset['TF'].unique()))
        n_tfs_common = len(TFs_common)
        n_tfs_common_list.append(n_tfs_common)
        
        df_tissue_subset = df_tissue_subset[df_tissue_subset['TF'].isin(TFs_common)]
        df_metanet_subset = df_metanet_subset[df_metanet_subset['TF'].isin(TFs_common)]

        print("> saving panda subset network")
        panda_tissue_subset_output_path = os.path.join(args.metanet_path, tissue_metanet, "benchmark_networks", "panda_subset")
        os.makedirs(panda_tissue_subset_output_path, exist_ok=True)
        df_tissue_subset.to_csv(os.path.join(panda_tissue_subset_output_path, f"{tissue_metanet}_panda_subset_stableID.tsv"), sep="\t", index=False, header=False)
        print("> saving metanet panda subset network")
        metanet_subset_panda_output_path = os.path.join(args.metanet_path, tissue_metanet, "METANet", "xgboost_subset_panda")
        os.makedirs(metanet_subset_panda_output_path, exist_ok=True)
        df_metanet_subset.to_csv(os.path.join(metanet_subset_panda_output_path, f"{tissue_metanet}_xgboost_subset_panda_stableID.tsv"), sep="\t", index=False, header=False)

        TFs_common_specific = set(df_metanet['TF'].unique()).intersection(set(df_tissue_specific['TF'].unique()))
        df_tissue_specific = df_tissue_specific[df_tissue_specific['TF'].isin(TFs_common_specific)]
        
        print("> saving tissue-specific panda subset network")
        panda_tissue_specific_subset_output_path = os.path.join(args.metanet_path, tissue_metanet, "benchmark_networks", "panda_tissue_specific")
        os.makedirs(panda_tissue_specific_subset_output_path, exist_ok=True)
        df_tissue_specific.to_csv(os.path.join(panda_tissue_specific_subset_output_path, f"{tissue_metanet}_panda_tissue_specific_stableID.tsv"), sep="\t", index=False, header=False)
        
        
        # calculate max targets per TF
        max_targets_per_tf = int(len(df_tissue_subset) / df_tissue_subset['TF'].nunique())
        print(f"\t> max {max_targets_per_tf} targets per TF")
    except FileNotFoundError:
        print(f"\t>{tissue_panda} not in METANets")


print("DONE")
