import pandas as pd
import numpy as np
import os
 
promoter_scores_df = pd.read_csv('/scratch/mblab/d.p.ruskin/human/MEME/output/sum_hits_by_chr/merged_promoter_scores.txt', sep="\t", comment="#")
maximum = promoter_scores_df['promoter_score'].max()
minimum = promoter_scores_df['promoter_score'].min()
def normalize(score):
    return (score - minimum)/(maximum - minimum)
promoter_scores_df['promoter_score'] = promoter_scores_df['promoter_score'].apply(normalize)

promoter_scores_df.to_csv('/scratch/mblab/d.p.ruskin/human/MEME/output/sum_hits_by_chr/merged_normalized_promoter_scores.txt', sep="\t", index=False)
