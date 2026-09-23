"""G1 validation: WHO 2023 grade-1/2 rpoB RIF entries vs locked 8.0 A rule; G1b on grade 4/5."""
import re, pandas as pd
AA3 = 'Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split()
d = pd.read_csv('results/residue_distances.tsv', sep='\t').dropna(subset=['who_pos'])
dist = dict(zip(d.who_pos.astype(int), d.min_dist))
w = pd.read_csv('data/who2023_rpoB_RIF.tsv', sep='\t')
w['grade'] = w['FINAL CONFIDENCE GRADING'].str[0].astype(int)
def span(v):
    m = v.split('p.', 1)[1] if 'p.' in v else ''
    pos = [int(x) for x in re.findall(r'[A-Z][a-z]{2}(\d+)', m)]
    if not pos: return None
    if 'ins' in m and 'delins' not in m and len(pos) == 2: return pos  # flanks
    return list(range(min(pos), max(pos) + 1))
def score(v):
    s = span(v)
    if s is None: return None, None, None
    ds = [dist[p] for p in s if p in dist]
    return (min(ds) if ds else float('nan')), min(s), max(s)
w[['score', 'pos_start', 'pos_end']] = w.variant.apply(lambda v: pd.Series(score(v)))
w['protein_level'] = w.effect.isin(['missense_variant', 'inframe_deletion', 'inframe_insertion', 'stop_gained', 'frameshift'])
def cls(x):
    if pd.isna(x): return 'unresolved/unscorable'
    return 'contact' if x <= 4.5 else ('shell' if x <= 8.0 else 'distal')
w['class'] = w.score.apply(cls)
w['plausible'] = w.score.le(8.0)
w.to_csv('results/who_rpoB_scored.tsv', sep='\t', index=False)
g12 = w[w.grade <= 2]
print('G1 denominator (all grade 1/2 rows):', len(g12), g12.effect.value_counts().to_dict())
pl = g12[g12.protein_level]
print('G1 protein-level entries:', len(pl), 'plausible:', pl.plausible.sum(), 'rate %.4f' % pl.plausible.mean())
print(pl['class'].value_counts().to_dict())
print('non-protein-level grade1/2:', g12[~g12.protein_level][['variant', 'effect']].to_string())
print('discrepancies:'); print(pl[~pl.plausible][['variant', 'grade', 'score', 'class']].to_string())
for lab, sub in [('grade3', w[(w.grade == 3) & w.protein_level]), ('grade4/5', w[(w.grade >= 4) & w.protein_level])]:
    print(lab, len(sub), 'plausible rate %.4f' % sub.plausible.mean(), sub['class'].value_counts().to_dict())
