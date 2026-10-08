#!/usr/bin/env python3
"""02_validate.py - apply the locked rule to WHO 2023 inhA entries; G1, G1b, discrepancies,
promoter line. All rule parameters are fixed in GATES_LOCKED.md + GATES_ADDENDUM_1.md."""
import re, os, pandas as pd
D = os.path.join(os.path.dirname(__file__), '..', 'data')
R = os.path.join(os.path.dirname(__file__), '..', 'results')
AA3 = 'Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split()
M3 = {a[0].upper()+a[1:].lower(): a for a in AA3}  # 3-letter -> itself
S3 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H',
      'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
CHARGE = {'R':1,'K':1,'H':0,'D':-1,'E':-1}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
ANCHORS = {94, 158, 165}

d = pd.read_csv(os.path.join(R, 'residue_distances.tsv'), sep='\t')
dist = dict(zip(d.who_pos.astype(int), d.min_dist))
rsa = dict(zip(d.who_pos.astype(int), d.rsa))

w = pd.read_csv(os.path.join(D, 'who2023_inhA_fabG1_INH_ETH.tsv'), sep='\t', dtype=str)
w['grade'] = w['FINAL CONFIDENCE GRADING'].str[0].astype(int)

def parse_protein(v):
    """Return (span list, wt1, mut1) from inhA_p.HGVS, or (None,None,None)."""
    if 'p.' not in v: return None, None, None
    m = v.split('p.', 1)[1]
    mm = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', m)
    if mm: return [int(mm.group(2))], S3.get(mm.group(1)), S3.get(mm.group(3))
    pos = [int(x) for x in re.findall(r'[A-Z][a-z]{2}(\d+)', m)]
    if not pos: return None, None, None
    if 'ins' in m and 'delins' not in m and len(pos) == 2: return pos, None, None
    return list(range(min(pos), max(pos) + 1)), None, None

w[['span','wt','mut']] = w.variant.apply(lambda v: pd.Series(parse_protein(v)))
w['pos_start'] = w.span.apply(lambda s: min(s) if s else None)
w['pos_end'] = w.span.apply(lambda s: max(s) if s else None)

def score_span(s):
    if not s: return float('nan')
    ds = [dist[p] for p in s if p in dist]
    return min(ds) if ds else float('nan')
w['score'] = w.span.apply(score_span)

protein_level = {'missense_variant','stop_gained','frameshift','inframe_deletion','inframe_insertion'}
w['protein_level'] = w.effect.isin(protein_level)

def classify(row):
    """Locked rule R1-R3 -> STRUCT-R / STRUCT-neutral (R0 has no instances in this subset)."""
    if not row.protein_level or pd.isna(row.score): return None
    pos = row.pos_start
    r1 = row.score <= 8.0
    r2 = any(p in ANCHORS for p in row.span)
    r3 = False
    if row.wt and row.mut and row.wt in VOL and row.mut in VOL:
        buried = rsa.get(pos, 1.0) < 0.20
        destab = (CHARGE.get(row.wt,0) != CHARGE.get(row.mut,0)) or \
                 ('P' in (row.wt, row.mut) and row.wt != row.mut) or \
                 abs(VOL[row.mut]-VOL[row.wt]) > 80
        r3 = buried and destab
    return {'R1': r1, 'R2': r2, 'R3': r3, 'STRUCT_R': r1 or r2 or r3}
cls = pd.DataFrame([classify(row) if classify(row) is not None else {'R1': None, 'R2': None, 'R3': None, 'STRUCT_R': None} for row in w.itertuples()], index=w.index)
w = pd.concat([w, cls], axis=1)
w['class'] = w.score.apply(lambda x: 'unresolved' if pd.isna(x) else ('contact' if x <= 4.5 else ('shell' if x <= 8.0 else 'distal')))
w.to_csv(os.path.join(R, 'who_inhA_scored.tsv'), sep='\t', index=False)

log = []
P = lambda *a: (log.append(' '.join(str(x) for x in a)), print(*a))
P('=== inhA G1 validation vs WHO 2023 (locked rule) ===')
P('WHO subset rows:', len(w), 'grades:', w.grade.value_counts().sort_index().to_dict())
allg12 = w[w.grade <= 2]
P('all grade-1/2 rows:', len(allg12), 'by effect:', allg12.effect.value_counts().to_dict())
for _, r in allg12.iterrows():
    P('  ', r.variant, '|', r.drug, '| grade', r.grade, '|', r.effect)

g12 = w[(w.grade <= 2) & w.protein_level].copy()
# union across drugs: unique variants (a variant may be grade 1/2 for INH and/or ETH)
g12u = g12.groupby('variant').agg(drug=('drug', lambda s: '+'.join(sorted(set(s)))),
                                  grade=('grade','min'), score=('score','first'),
                                  STRUCT_R=('STRUCT_R','first'), span=('span','first'),
                                  clas=('class','first'), R1=('R1','first'), R2=('R2','first'), R3=('R3','first')).reset_index()
P()
P('G1 denominator (unique grade-1/2 protein-altering variants, union across drugs):', len(g12u))
for _, r in g12u.iterrows():
    P(' ', r.variant, '|', r.drug, '| grade', r.grade, '| dist', r.score, '|', r.clas,
      '| R1', r.R1, 'R2', r.R2, 'R3', r.R3, '->', 'STRUCT-R' if r.STRUCT_R else 'STRUCT-NEUTRAL')
hits = int((g12u.STRUCT_R == True).sum())
P()
P('G1 RESULT: %d/%d = %.1f%% STRUCT-R  (gate >= 90%%)' % (hits, len(g12u), 100*hits/len(g12u)))
P('G1:', 'PASS' if hits/len(g12u) >= 0.90 else 'FAIL')
P()
per = g12.groupby('drug').apply(lambda s: ((s.STRUCT_R == True).sum(), len(s)), include_groups=False)
P('per-drug (row-level):', {k: f'{a}/{b}' for k,(a,b) in per.items()})
P()
disc = g12u[g12u.STRUCT_R != True]
P('DISCREPANCIES (grade-1/2 labeled STRUCT-neutral):', len(disc))
for _, r in disc.iterrows():
    P(' ', r.variant, r.drug, 'grade', r.grade, 'dist', r.score, r.clas)
P()
g45 = w[(w.grade >= 4) & w.protein_level]
P('G1b specificity: grade-4/5 protein-altering rows:', len(g45),
  'STRUCT-R rate %.1f%%' % (100*(g45.STRUCT_R == True).mean()))
g3 = w[(w.grade == 3) & w.protein_level]
P('grade-3 (uncertain, reported separately):', len(g3), 'rows, STRUCT-R rate %.1f%%' % (100*g3.STRUCT_R.mean()))
P()
up = w[w.effect == 'upstream_gene_variant'].copy()
up['upos'] = up.mutation.str.extract(r'c\.-(\d+)').astype(float)
prom = up[up.upos <= 50]
P('PROMOTER LINE (descriptive, not gated): upstream entries in -1..-50 class:', len(prom))
P('  grade distribution:', prom.grade.value_counts().sort_index().to_dict())
for _, r in prom[prom.grade <= 2].iterrows():
    P('  grade-1/2 promoter entry:', r.variant, '|', r.drug, '| grade', r.grade)
syn = w[w.effect == 'synonymous_variant']
P()
P('synonymous (context):', len(syn), 'rows, grades:', syn.grade.value_counts().sort_index().to_dict())
open(os.path.join(R, 'G1_validation_log.txt'), 'w').write('\n'.join(x if isinstance(x,str) else str(x) for x in log))
