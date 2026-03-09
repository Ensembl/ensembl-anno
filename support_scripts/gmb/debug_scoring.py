import pyranges as pr
import pandas as pd
from config import load_config
from scoring import select_isoforms, score_model, merge_identical_models

config = load_config('configs/fungi_default.yaml')

df = pr.read_gff3('z_tritici/helixer_remapped.gff3').df
e_mask = df['Feature'].isin(['exon', 'CDS'])
df.loc[e_mask, 'transcript_id'] = df.loc[e_mask, 'Parent']
df['Source'] = 'Helixer'

cluster = df[df['Feature'] == 'exon'].head(50)

models = []
for tid, tdf in cluster.groupby('transcript_id'):
    s = score_model(tdf, config.scoring)
    if s is not None:
        models.append(s)

print(f"Generated {len(models)} models.")
for m in models:
    print(f"Model ID: {m['id']} Score: {m['score']} Exons: {m['rep']['exon_count']}")

merged = merge_identical_models(models)
print(f"Merged into {len(merged)} models.")

candidates = []
for s in merged.values():
    keep = False
    if 'Helixer' in s['sources'] and config.scoring.keep_helixer_without_support:
        keep = True
    print(f"Candidate {s['id']} kept? {keep} Sources: {s['sources']}")
