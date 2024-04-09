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



## JASPAR 2022 Motifs
see `data/JASPAR2022`

[https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz)
> Castro-Mondragon, J. et al., JASPAR 2022: the 9th release of the open-access
> database of transcription factor binding profiles. Nucleic Acids Research,
> Volume 50, Issue D1, 7 January 2022, Pages D165–D173,
> https://doi.org/10.1093/nar/gkab1113

The JASPAR 2022 motifs were downloaded from the MEME Motif databases which were
already converted to the MEME format required by FIMO. Human TF motifs were
subsetted for from the JASPAR (NON-REDUNDANT) CORE (2022) vertebrates database.
Dimers were removed.


## Enhancer Atlas 2.0
see `data/enhanceratlas2.0`

[http://www.enhanceratlas.org/data/download/species_enh_bed.tar.gz](http://www.enhanceratlas.org/data/download/species_enh_bed.tar.gz)
> Gao, T. et al, EnhancerAtlas: a resource for enhancer annotation and analysis
> in 105 human cell/tissue types. Bioinformatics 2016; doi:
> 10.1093/bioinformatics/btw495.

I downloaded the dataset called "Enhancers of all species by bed format" under
Download/Download enhancers. Human enhancers were selected and corresponding
enhancer-gene interactions were downloaded individually from Download/Download
enhancer-gene interactions.

## ENCODE TF Perturbations
see `data/ENCODE`

[https://www.encodeproject.org/](https://www.encodeproject.org/)
> The ENCODE Project Consortium. An integrated encyclopedia of DNA elements in
> the human genome. Nature 489, 57–74 (2012).
> https://doi.org/10.1038/nature11247

> Luo Y., et al,  New developments on the Encyclopedia of DNA Elements (ENCODE)
> data portal.  Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020,
> Pages D882–D889, https://doi.org/10.1093/nar/gkz1062

> Hitz, B. C., et al,  The ENCODE Uniform Analysis Pipelines. bioRxiv
> 2023.04.04.535623; doi: https://doi.org/10.1101/2023.04.04.535623

I downloaded the accessions in `see/ENCODE/file_download_metadata_lamberts.tsv`
which contain accessions that meet the criteria: CRISPRi RNA-seq, shRNA RNA-seq,
siRNA RNA-seq with replicates for K562 or HepG2 homo sapiens biosamples.

## Dealing with Motif Overlaps
see `scripts/motif_overlaps`

## Connecting TF Motifs to Promoters
see `scripts/motif_to_promoter`

## Creating TF to TG Matrices
see `scripts/tf_to_tg`
