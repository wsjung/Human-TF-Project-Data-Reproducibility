# Human-TF-Project-Data-Reproducibility 
This project stores the data, scripts, citation and links which describe the source and transformations on externally
produced data which is added to the Human TF Network project 

## FANTOM5 CAGE
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

## Marbach TF networks
see `data/Marbach`
> Marbach D, Lamparter D, Quon G, Kellis M, Kutalik Z, Bergmann S:
> Tissue-specific regulatory circuits reveal variable modular perturbations
> across complex diseases. Nat Methods 2016, 13(4):366-370.
> https://doi.org/10.1038/nmeth.3799

## PANDA TF networks
see `data/PANDA`
> Sonawane AR, Platig J, Fagny M, Chen C-Y, Paulson JN, Lopes-Ramos CM, DeMeo
> DL, Quackenbush J, Glass K, Kuijjer ML: Understanding Tissue-Specific Gene
> Regulation. Cell Reports 2017, 21(4):1077-1088.
> https://doi.org/10.1016/j.celrep.2017.10.001

## Dealing with Motif Overlaps
see `scripts/motif_overlaps`

## Connecting TF Motifs to Promoters
see `scripts/motif_to_promoter`

## Creating TF to TG Matrices
For TF->TG matrices based on FANTOM5, see `scripts/tf_to_tg`.

For TF->TG matrices based on EnhancerAtlas, see `scripts/tf_to_tg_matrix_enhanceratlas`.

## Creating Binding labels based on REMAP 2020 peaks
see `scripts/remap_labels`

## NP3 evaluation scripts
see `scripts/evaluation/` for evaluation scripts of the NP3 network and its input feature networks.

Note: NP3 evaluation scripts are based on
[NET-evaluation](https://github.com/BrentLab/NET-evaluation).

## Running BART on tissue-specific data
see `scripts/bart`

## Running BART on concatenated (tissue-agnostic) data
see `scripts/bart_concat`

## Running GENIE3
see `scripts/genie3`

## Running WGCNA
see `scripts/wgcna`

## Example eQTL evaluation script (36 human networks and tissue-specific eQTL)
see `scripts/evaluation/eQTL_analysis_original_filtered`

## Example eQTL evaluation script (1 aggregate network and tissue-specific eQTL)
see `scripts/evaluation/eQTL_analysis_aggregate_filtered`
