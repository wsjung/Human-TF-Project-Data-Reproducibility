# REMAP 2020 binding score processing
=====================================

Scripts for computing binding scores for each TF-TG edge.
* Note that the actual calculation of binary labels for XGBoost classification
  modeling can be found in `../metanet_modeling/xgboost_5.py`.

## Steps
--------

1. Filtering of significant (FDR q<=0.01) binding peaks
    * `filter_sig_qval_peaks_1.py`

2. Mapping of binding peaks to FANTOM5 regulatory elements
    * `map_remap_to_fantom5_elements_2.py`
    * Note that we only keep regulatory elements that overlap with remap peaks
      if the ovelapping region is at least half the length of the remap peak.

3. Takes the maximum peak score for overlapping remap peaks and calculates a
   score for each promoter/enhancer
    * `max_remap_peak_sum_scores_3.py`

4. Annotates promoters and enhancers with target genes
    * `annotate_fantom5_remap_peaks_with_genes_4.py`
    * Annotations are based on FANTOM5.

5. Calculate TF-TG binding score by summing all CRE paths
    * `sum_tf_tg_scores_5.py`
    * A score for each gene is derived by summing the scores of all regulatory
      elements annotated for the gene.
