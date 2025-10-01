# METANet XGBoost modeling
==========================

Scripts for preparing the input files for training XGBoost models and making
predictions to infer METANets.

## Input
--------

These scripts assume previous data processing steps have been followed and the
output files for the tissue-specific and tissue-aggregate data are saved in the
following directory structure:
```
raw_data
├── adipose_subcutaneous
│   ├── LASSO.txt
│   ├── MARBACH.tsv
│   └── net_bart.tsv
├── adipose_visceral_omentum
│      .
│      .
│      .
└── whole_blood
    ├── LASSO.txt
    ├── MARBACH.tsv
    └── net_bart.tsv
raw_data_shared
├── BART_all_gtex_tissues_wide.tsv
├── LASSO_all_gtex_tissues_wide.tsv
└── Remap_binding_scores.bed
```
* `raw_data` contains tissue-specific data.
* `raw_data_shared` contains tissue-aggregate data.


## Steps
--------

1. Find common TFs and genes and binarize binding data to positive / negative labels.
    * `process_common_tfs_and_genes_1.py`
    * Creates an output file listing all feature values and labels for each TF-gene pair

2. Feature transformation.
    * `feature_transformation_2.R`
    * This script performs the -log rank transformation + min-max normalization
      for  each feature.

3. Extract features to individual files.
    * `extract_features_to_files_3.py`
    * This script extracts the feature and label columns from the above
      step's output file to individual files.

4. Create cross-validation (CV) folds.
    * `create_cv_folds_4.py`
    * This script constructs the (k=10) cross-validation folds, stratified by
      binding label proportions.

5. Train XGBoost models and make predictions.
    * `xgboost_5.py`
    * This script carries out the training of the XGBoost models, one for each
      of the 10 cross-validation folds. Each trained XGBoost model is saved as
      a JSON file and used to make predictions on the test fold data.


## Output
---------
This pipeline outputs
1. `CV_folds` folder which stores the cross-validation data from
   `create_cv_folds_4.py`.
2. `xgboost/CV` folder which stores the model JSON file, raw predictions, tuned
   hyperparameters, etc.
3. `xgboost/${tissue}_xgboost_full.tsv` and `xgboost/${tissue}_xgboost.tsv`
   which are the predicted probabilities for all TF-gene edges. These are the
   METANets.

```
results
├── adipose_subcutaneous
│   ├── 10.0_threshold_metanet
│   │   ├── CV_folds
│   │   │   ├── fold0_test_data.txt
│   │   │   ├── fold0_train_data.txt
│   │   │   ├── fold1_test_data.txt
│   │   │   ├── fold1_train_data.txt
│   │   │   ├──        ...
│   │   │   ├── fold9_test_data.txt
│   │   │   └── fold9_train_data.txt
│   │   └── xgboost
│   │       ├── adipose_subcutaneous_xgboost_full.tsv
│   │       ├── adipose_subcutaneous_xgboost.tsv
│   │       ├── CV
│   │       │   ├── fold0_auprc.png
│   │       │   ├── fold0_auroc.png
│   │       │   ├── fold0_hyperparams.tsv
│   │       │   ├── fold0_model.json
│   │       │   ├── fold0_test_preds.tsv
│   │       │   ├──       ...
│   │       │   ├── fold9_auprc.png
│   │       │   ├── fold9_auroc.png
│   │       │   ├── fold9_hyperparams.tsv
│   │       │   ├── fold9_model.json
│   │       │   ├── fold9_test_preds.tsv
```

