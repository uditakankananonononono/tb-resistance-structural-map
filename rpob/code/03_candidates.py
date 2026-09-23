"""G2: CRyPTIC solo-isolate MIC evidence for rpoB variants; novel candidates = grade 3 / uncatalogued,
<=8.0 A from RFP, >=3 isolates carrying the variant as their only rpoB protein-altering variant, median MIC > ECOFF.
ECOFF = 0.5 mg/L (highest MIC labelled S in CRyPTIC BINARY_PHENOTYPE; R starts at 1.0)."""
import re, numpy as np, pandas as pd
A3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','Z':'Ter','!':'Ter'}
ECOFF = 0.5
dist = pd.read_csv('results/residue_distances.tsv', sep='\t').dropna(subset=['who_pos'])
dist = dict(zip(dist.who_pos.astype(int), dist.min_dist))
who = pd.read_csv('results/who_rpoB_scored.tsv', sep='\t')
gmap = dict(zip(who.variant.str.replace('rpoB_', '', regex=False), who.grade))
e = pd.read_csv('data/cryptic_rpoB_RIF_effects.tsv.gz', sep='\t')
def parse(m):
    x = re.fullmatch(r'([A-Z!])(-?\d+)([A-Z!])', m)
    if x:
        a, p, b = x.group(1), int(x.group(2)), x.group(3)
        if a == b or 'X' in (a, b): return None  # synonymous or null call
        return dict(kind='missense', codon_start=p, codon_end=p, hgvs=f'p.{A3[a]}{p}{A3[b]}')
    x = re.fullmatch(r'(-?\d+)_(ins|del)_([a-z]+)', m)
    if x:
        p, t, s = int(x.group(1)), x.group(2), x.group(3)
        if p < 1: return None  # upstream
        c0 = (p - 1) // 3 + 1; c1 = (p + len(s) - 1) // 3 + 1 if t == 'del' else c0 + 1
        return dict(kind=('inframe_' if len(s) % 3 == 0 else 'frameshift_') + t, codon_start=c0, codon_end=c1, hgvs=None)
    return None  # promoter/nucleotide calls
rows = []
for r in e.itertuples():
    q = parse(r.MUTATION)
    if q: rows.append(dict(UNIQUEID=r.UNIQUEID, MUTATION=r.MUTATION, **q))
v = pd.DataFrame(rows)
ev = e.set_index(['UNIQUEID', 'MUTATION']).EVIDENCE
def who_hgvs(u, m, h):
    s = str(ev.get((u, m), ''))
    x = re.search(r"'WHO HGVS': 'rpoB_(p\.[^']+)'", s)
    return x.group(1) if x else h
v['hgvs'] = [who_hgvs(a, b, c) for a, b, c in zip(v.UNIQUEID, v.MUTATION, v.hgvs)]
v['who_grade'] = v.hgvs.map(gmap)
v['score'] = [min([dist[p] for p in range(a, b + 1) if p in dist] or [np.nan]) for a, b in zip(v.codon_start, v.codon_end)]
n_per = v.groupby('UNIQUEID').size()
v['solo'] = v.UNIQUEID.map(n_per).eq(1)
p = pd.read_csv('data/cryptic_RIF_phenotypes.tsv.gz', sep='\t')
def mic(s):
    s = str(s)
    if s.startswith('<='): return float(s[2:]) / 2
    if s.startswith('>'): return float(s[1:]) * 2
    return float(s)
p['mic'] = p.MIC.map(mic)
p = p[['UNIQUEID', 'mic', 'BINARY_PHENOTYPE', 'PHENOTYPE_QUALITY']]
vs = v[v.solo].merge(p, on='UNIQUEID')
# wild type: phenotyped samples in EFFECTS with no rpoB protein-altering variant
wt_ids = set(e.UNIQUEID) - set(v.UNIQUEID)
wt = p[p.UNIQUEID.isin(wt_ids)]
key = vs.MUTATION
agg = vs.groupby('MUTATION').agg(hgvs=('hgvs', 'first'), kind=('kind', 'first'), codon=('codon_start', 'first'),
    who_grade=('who_grade', 'first'), score=('score', 'first'), n_solo=('UNIQUEID', 'nunique'),
    median_mic=('mic', 'median'), frac_R=('BINARY_PHENOTYPE', lambda s: (s == 'R').mean()),
    n_high=('PHENOTYPE_QUALITY', lambda s: (s == 'HIGH').sum())).reset_index()
agg['class'] = pd.cut(agg.score, [-1, 4.5, 8.0, 1e9], labels=['contact', 'shell', 'distal']).astype(str)
agg.loc[agg.score.isna(), 'class'] = 'unresolved'
agg['status'] = agg.who_grade.map({1: 'grade1', 2: 'grade2', 3: 'grade3', 4: 'grade4', 5: 'grade5'}).fillna('uncatalogued')
agg = agg.sort_values(['n_solo'], ascending=False)
agg.to_csv('results/cryptic_solo_variant_mic.tsv', sep='\t', index=False)
novel = agg[agg.status.isin(['grade3', 'uncatalogued']) & (agg.score <= 8.0)]
cand = novel[(novel.n_solo >= 3) & (novel.median_mic > ECOFF)]
thin = novel[(novel.n_solo < 3) & (novel.median_mic > ECOFF)]
cand.to_csv('results/G2_candidates.tsv', sep='\t', index=False); thin.to_csv('results/G2_thin.tsv', sep='\t', index=False)
novel.to_csv('results/G2_all_plausible_novel.tsv', sep='\t', index=False)
print('solo phenotyped isolates', vs.UNIQUEID.nunique(), 'WT phenotyped', len(wt), 'WT frac R %.4f' % (wt.BINARY_PHENOTYPE == 'R').mean(), 'WT median MIC', wt.mic.median())
print('variants with solo MIC:', len(agg)); print(agg.groupby('status').size().to_dict())
print('plausible novel variants with MIC:', len(novel), '| G2 candidates:', len(cand), '| thin:', len(thin))
print(cand[['MUTATION', 'hgvs', 'status', 'score', 'class', 'n_solo', 'n_high', 'median_mic', 'frac_R']].to_string(index=False))
print('--- plausible novel below ECOFF (n>=3):'); print(novel[(novel.n_solo >= 3) & (novel.median_mic <= ECOFF)][['MUTATION', 'status', 'score', 'n_solo', 'median_mic', 'frac_R']].to_string(index=False))
print('--- distal novel above ECOFF n>=3 (outside rule):'); print(agg[agg.status.isin(['grade3', 'uncatalogued']) & (agg.score > 8) & (agg.n_solo >= 3) & (agg.median_mic > ECOFF)][['MUTATION', 'status', 'score', 'n_solo', 'median_mic', 'frac_R']].to_string(index=False))
