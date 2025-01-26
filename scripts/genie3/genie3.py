import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import argparse

# Function to compute feature importances
def compute_feature_importances(estimator):
    importances = [e.tree_.compute_feature_importances(normalize=False) for e in estimator.estimators_]
    importances = np.array(importances)
    return np.sum(importances, axis=0) / len(estimator)

# Function to run GENIE3 for a single gene
def GENIE3_single(expr_data, output_idx, input_idx, tree_method, K, ntrees):
    ngenes = expr_data.shape[1]
    output = expr_data[:, output_idx]
    output = output / np.std(output)
    input_idx = input_idx[:]
    if output_idx in input_idx:
        input_idx.remove(output_idx)
    expr_data_input = expr_data[:, input_idx]

    if (K == 'all') or (isinstance(K, int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K

    if tree_method == 'RF':
        #treeEstimator = RandomForestRegressor(n_estimators=ntrees, max_features=max_features)
        treeEstimator = RandomForestRegressor(n_estimators=ntrees, max_features=max_features, random_state=0) # Added by Daniel to get same results each time

    treeEstimator.fit(expr_data_input, output)
    feature_importances = compute_feature_importances(treeEstimator)
    vi = np.zeros(ngenes)
    vi[input_idx] = feature_importances
    return vi

# Function to run GENIE3 for all genes
def GENIE3(expr_data, gene_names=None, regulators='all', tree_method='RF', K='sqrt', ntrees=1000, nthreads=1):
    ngenes = expr_data.shape[1]
    if gene_names is not None and len(gene_names) != ngenes:
        raise ValueError('gene_names must be a list of length equal to the number of columns in expr_data')
    
    input_idx = list(range(ngenes)) if regulators == 'all' else [i for i, gene in enumerate(gene_names) if gene in regulators]
    
    VIM = np.zeros((ngenes, ngenes))
    for i in range(ngenes):
        vi = GENIE3_single(expr_data, i, input_idx, tree_method, K, ntrees)
        VIM[i, :] = vi
    
    return VIM.T

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run GENIE3 for a batch of genes.')
parser.add_argument('--input_file', required=True, help='Path to the input CSV file containing expression data.')
parser.add_argument('--batch_file', required=True, help='Path to the batch file containing gene names.')
parser.add_argument('--output_dir', required=True, help='Directory to save the output CSV files.')
args = parser.parse_args()

# Load the input data
data = pd.read_csv(args.input_file)

# Extract gene names (column headers)
gene_names = data.columns.tolist()

# Get the expression matrix
expr_data = data.values

# Read gene names from the batch file
with open(args.batch_file, 'r') as f:
    batch_genes = [line.strip() for line in f]

# Run GENIE3 for each gene in the batch
for gene in batch_genes:
    if gene in gene_names:
        output_idx = gene_names.index(gene)
        feature_importances = GENIE3_single(expr_data, output_idx, list(range(expr_data.shape[1])), 'RF', 'sqrt', 1000)
        
        # Create a DataFrame with gene IDs and their feature importances
        fi_df = pd.DataFrame({
            'gene_id': gene_names,
            'feature_importance': feature_importances
        })
        
        # Save the DataFrame to a CSV file named after the gene_id
        fi_df.to_csv(f"{args.output_dir}/{gene}.csv", index=False)

print(f"Feature importance calculation and saving completed for batch: {args.batch_file}.")

