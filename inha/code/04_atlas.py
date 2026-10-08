#!/usr/bin/env python3
"""04_atlas.py - per-residue InhA structural atlas: distances, RSA, WHO grade, CRyPTIC evidence."""
import os, pandas as pd, numpy as np
R = os.path.join(os.path.dirname(__file__), '..', 'results')
d = pd.read_csv(os.path.join(R, 'residue_distances.tsv'), sep='\t')
who = pd.read_csv(os.path.join(R, 'who_inhA_scored.tsv'), sep='\t')
w = who[who.protein_level & who.pos_start.notna()].copy()
w['pos_start'] = w.pos_start.astype(int)
w = w.rename(columns={'pos_start':'who_pos'})
var_by_pos = w.groupby('who_pos').agg(who_variants=('variant', lambda s: ';'.join(sorted(set(s)))),
    best_grade=('grade','min'), drugs=('drug', lambda s: '+'.join(sorted(set(s))))).reset_index()
at = d.merge(var_by_pos, on='who_pos', how='left')
at['mech_class'] = pd.cut(at.min_dist, [-1,4.5,8.0,1e9], labels=['direct contact','pocket shell','distal']).astype(str)
for drug in ('INH','ETH'):
    m = pd.read_csv(os.path.join(R, f'cryptic_solo_variant_mic_{drug}.tsv'), sep='\t')
    m = m.drop_duplicates(subset=['codon']).set_index('codon')
    at[f'n_solo_{drug}'] = at.who_pos.map(m.n_solo)
    at[f'median_mic_{drug}'] = at.who_pos.map(m.median_mic)
    at[f'fracR_{drug}'] = at.who_pos.map(m.frac_R)
at = at.sort_values('min_dist')
at.to_csv(os.path.join(R, 'inhA_structural_atlas.tsv'), sep='\t', index=False)
print('atlas rows:', len(at))
print(at[at.best_grade.notna()][['who_pos','resname','min_dist','mech_class','best_grade','who_variants','drugs']].to_string(index=False))
