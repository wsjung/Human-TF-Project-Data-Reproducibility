# Human-TF-Project-Data-Reproducibility 
This project stores the data, scripts, citation and links which describe the source and transformations on externally
produced data which is added to the Human TF Network project 

## Data

### FANTOM5 CAGE
see `data/FANTOM5`

[https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/](https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/)
> The FANTOM Consortium and the RIKEN PMI and CLST (DGT). A promoter-level
> mammalian expression atlas. Nature 507, 462–470 (2014).
> https://doi.org/10.1038/nature13182

> Andersson, R., Gebhard, C., Miguel-Escalada, I. et al. An atlas of active
> enhancers across human cell types and tissues. Nature 507, 455–461 (2014).
> https://doi.org/10.1038/nature12787

> Erik Arner et al., Transcribed enhancers lead waves of coordinated
> transcription in transitioning mammalian
> cells.Science347,1010-1014(2015).DOI:10.1126/science.1259418 

### Marbach TF networks
see `data/Marbach`
> Marbach D, Lamparter D, Quon G, Kellis M, Kutalik Z, Bergmann S:
> Tissue-specific regulatory circuits reveal variable modular perturbations
> across complex diseases. Nat Methods 2016, 13(4):366-370.
> https://doi.org/10.1038/nmeth.3799

### PANDA TF networks
see `data/PANDA`
> Sonawane AR, Platig J, Fagny M, Chen C-Y, Paulson JN, Lopes-Ramos CM, DeMeo
> DL, Quackenbush J, Glass K, Kuijjer ML: Understanding Tissue-Specific Gene
> Regulation. Cell Reports 2017, 21(4):1077-1088.
> https://doi.org/10.1016/j.celrep.2017.10.001

## Scripts

### Processing REMAP 2020 TF ChIP-seq binding scores
see `scripts/remap_labels`

### Expression filtering
see `scripts/expression_filtering`

### LASSO expression feature construction
see `scripts/lasso`

### BART expression feature construction
see `scripts/bart`

### METANet XGBoost modeling
> Scripts for preprocessing, modeling, and predictions using XGBoost.
see `metanet_modeling`

### Subsetting of PANDA networks
see `scripts/PANDA_subsetting`

### Network quality evaluation scripts
see `scripts/evaluation/network_quality/` for network quality evaluation
scripts.
> GO and GO-directness evaluation were done using
> [NET-evaluation](https://github.com/BrentLab/NET-evaluation).

### eQTL evaluation script
see `scripts/evaluation/tissue_specificity`

### Interpreting evaluation results
see `scripts/notebooks` for jupyter notebooks of interpreting evaluation
results.

### FISHNET application
see `scripts/fishnet`
