#!/usr/bin/env python3
"""03_candidates.py - G2: CRyPTIC solo-isolate MIC evidence for inhA variants.
Candidates = grade-3/uncatalogued coding variants labeled STRUCT-R by the locked rule,
sole inhA/fabG1 protein-altering or promoter variant in >=3 phenotyped isolates,
median MIC > drug ECOFF. ECOFF derived per drug: highest MIC labelled S (CRyPTIC binary)."""
import re, os, numpy as np, pandas as pd
D = os.path.join(os.path.dirname(__file__), '..', 'data')
R = os.path.join(os.path.dirname(__file__), '..', 'results')
A3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','Z':'Ter','!':'Ter'}
CHARGE = {'R':1,'K':1,'H':0,'D':-1,'E':-1}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
ANCHORS = {94, 158, 165}
dist_df = pd.read_csv(os.path.join(R, 'residue_distances.tsv'), sep='\t')
dist = dict(zip(dist_df.who_pos.astype(int), dist_df.min_dist))
rsa = dict(zip(dist_df.who_pos.astype(int), dist_df.rsa))
who = pd.read_csv(os.path.join(R, 'who_inhA_scored.tsv'), sep='\t')
gmap = dict(zip(who.variant.str.replace('inhA_', '', regex=False), who.grade))

e = pd.read_csv(os.path.join(D, 'cryptic_inhA_fabG1_INH_ETH_effects.tsv.gz'), sep='\t')
def parse(m):
    x = re.fullmatch(r'([A-Z!])(-?\d+)([A-Z!])', m)
    if x:
        a, p, b = x.group(1), int(x.group(2)), x.group(3)
        if a == b or 'X' in (a, b): return None
        if p < 1: return None
        return dict(kind='missense', codon_start=p, codon_end=p, hgvs=f'p.{A3[a]}{p}{A3[b]}', wt=a, mut=b)
    return None  # promoter / indel calls: excluded from coding candidate scoring
rows = []
for r in e.itertuples():
    q = parse(r.MUTATION)
    if q: rows.append(dict(UNIQUEID=r.UNIQUEID, MUTATION=r.MUTATION, GENE=r.GENE, **q))
v = pd.DataFrame(rows).drop_duplicates(subset=['UNIQUEID','MUTATION'])
ev = e.set_index(['UNIQUEID','MUTATION']).EVIDENCE
def who_hgvs(u, m, h):
    s = str(ev.get((u, m), ''))
    x = re.search(r"'WHO HGVS': 'inhA_(p\.[^']+)'", s)
    return x.group(1) if x else h
v['hgvs'] = [who_hgvs(a,b,c) for a,b,c in zip(v.UNIQUEID, v.MUTATION, v.hgvs)]
v['who_grade'] = v.hgvs.map(gmap)
v['score'] = [min([dist[p] for p in range(a, b+1) if p in dist] or [np.nan]) for a,b in zip(v.codon_start, v.codon_end)]
# full locked rule (R1-R3)
def rule(r):
    r1 = pd.notna(r.score) and r.score <= 8.0
    r2 = r.codon_start in ANCHORS
    r3 = False
    if r.wt in VOL and r.mut in VOL and pd.notna(r.score):
        buried = rsa.get(r.codon_start, 1.0) < 0.20
        destab = (CHARGE.get(r.wt,0) != CHARGE.get(r.mut,0)) or ('P' in (r.wt, r.mut) and r.wt != r.mut) or abs(VOL[r.mut]-VOL[r.wt]) > 80
        r3 = bool(buried and destab)
    return bool(r1 or r2 or r3)
v['STRUCT_R'] = v.apply(rule, axis=1)
# solo: only ONE inhA/fabG1 protein-altering variant in the isolate, and no upstream variants there
up_ids = set(e[e.MUTATION.str.match(r'^[a-z]-\d+[a-z]$', na=False)].UNIQUEID)
n_per = v.groupby('UNIQUEID').size()
v['solo'] = v.UNIQUEID.map(n_per).eq(1) & ~v.UNIQUEID.isin(up_ids)
p = pd.read_csv(os.path.join(D, 'cryptic_INH_ETH_phenotypes.tsv.gz'), sep='\t')
def mic(s):
    s = str(s)
    if s.startswith('<='): return float(s[2:]) / 2
    if s.startswith('>'): return float(s[1:]) * 2
    return float(s)
p['mic'] = p.MIC.map(mic)
out_all = {}
for drug in ('INH', 'ETH'):
    pd_ = p[p.DRUG == drug]
    s_mics = pd_.loc[pd_.BINARY_PHENOTYPE == 'S', 'mic']
    ecoff = s_mics.max()
    vs = v[v.solo].merge(pd_[['UNIQUEID','mic','BINARY_PHENOTYPE','PHENOTYPE_QUALITY']].drop_duplicates('UNIQUEID'), on='UNIQUEID')
    wt_ids = set(e.UNIQUEID) - set(v.UNIQUEID) - up_ids
    wt = pd_[pd_.UNIQUEID.isin(wt_ids)]
    agg = vs.groupby('MUTATION').agg(hgvs=('hgvs','first'), codon=('codon_start','first'),
        who_grade=('who_grade','first'), score=('score','first'), STRUCT_R=('STRUCT_R','first'),
        n_solo=('UNIQUEID','nunique'), median_mic=('mic','median'),
        frac_R=('BINARY_PHENOTYPE', lambda s: (s=='R').mean()),
        n_high=('PHENOTYPE_QUALITY', lambda s: (s=='HIGH').sum())).reset_index()
    agg['class'] = pd.cut(agg.score, [-1, 4.5, 8.0, 1e9], labels=['contact','shell','distal']).astype(str)
    agg.loc[agg.score.isna(), 'class'] = 'unresolved'
    agg['status'] = agg.who_grade.map({1:'grade1',2:'grade2',3:'grade3',4:'grade4',5:'grade5'}).fillna('uncatalogued')
    agg = agg.sort_values('n_solo', ascending=False)
    agg.to_csv(os.path.join(R, f'cryptic_solo_variant_mic_{drug}.tsv'), sep='\t', index=False)
    novel = agg[agg.status.isin(['grade3','uncatalogued']) & agg.STRUCT_R]
    cand = novel[(novel.n_solo >= 3) & (novel.median_mic > ecoff)]
    thin = novel[(novel.n_solo < 3) & (novel.median_mic > ecoff)]
    cand.to_csv(os.path.join(R, f'G2_candidates_{drug}.tsv'), sep='\t', index=False)
    thin.to_csv(os.path.join(R, f'G2_thin_{drug}.tsv'), sep='\t', index=False)
    novel.to_csv(os.path.join(R, f'G2_all_plausible_novel_{drug}.tsv'), sep='\t', index=False)
    print(f'== {drug} == ECOFF {ecoff} mg/L | solo phenotyped isolates {vs.UNIQUEID.nunique()} | WT n={len(wt)} fracR={(wt.BINARY_PHENOTYPE=="R").mean():.4f} medianMIC={wt.mic.median()}')
    print('  variants with solo MIC:', len(agg), agg.groupby('status').size().to_dict())
    print('  STRUCT-R novel:', len(novel), '| G2 candidates:', len(cand), '| thin:', len(thin))
    if len(cand): print(cand[['MUTATION','hgvs','status','score','class','n_solo','n_high','median_mic','frac_R']].to_string(index=False))
    if len(thin): print('  THIN:'); print(thin[['MUTATION','hgvs','status','score','class','n_solo','median_mic','frac_R']].to_string(index=False))
    out_all[drug] = dict(ecoff=ecoff, n_novel=len(novel), n_cand=len(cand), n_thin=len(thin),
                         wt_n=len(wt), wt_median=wt.mic.median(), wt_fracR=float((wt.BINARY_PHENOTYPE=='R').mean()))
print(out_all)
