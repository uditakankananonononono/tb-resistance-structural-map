"""Mechanism atlas: per-residue RIF contact details for every residue carrying a WHO grade-1/2 entry, plus exploratory distance-MIC analysis."""
import gemmi, numpy as np, pandas as pd
from scipy.stats import spearmanr
st = gemmi.read_structure('data/5UHC.cif.gz'); st.remove_hydrogens(); m = st[0]
rif = [(a.name, np.array([a.pos.x, a.pos.y, a.pos.z])) for ch in m for r in ch if r.name == 'RFP' for a in r]
# rifampicin moieties from the CCD RFP bond graph (data/RFP_ccd.cif): naphthofuranone core = C1-C14 + O1-O4, O12;
# hydrazone-piperazine tail = C43, N2, N3, N4, C38-C42; ansa bridge = C15-C37, N1, O5-O11
CORE = {f'C{i}' for i in range(1, 15)} | {'O1', 'O2', 'O3', 'O4', 'O12'}
TAIL = {'C43', 'N2', 'N3', 'N4'} | {f'C{i}' for i in range(38, 43)}
def moiety(n):
    return 'naphthofuranone core' if n in CORE else 'piperazine/hydrazone tail' if n in TAIL else 'ansa bridge'
d = pd.read_csv('results/residue_distances.tsv', sep='\t').dropna(subset=['who_pos'])
ch = m[d.chain.iloc[0]]
res = {r.seqid.num: r for r in ch}
who = pd.read_csv('results/who_rpoB_scored.tsv', sep='\t')
mic = pd.read_csv('results/cryptic_solo_variant_mic.tsv', sep='\t')
g12 = who[who.grade <= 2]
rows = []
for pos in sorted(set(int(x) for x in g12.pos_start.dropna()) | set(int(x) for x in g12.pos_end.dropna())):
    rr = res.get(pos + 6)
    if rr is None: continue
    best = (1e9, None, None)
    hb = []
    for a in rr:
        for n, x in rif:
            dd = np.linalg.norm(np.array([a.pos.x, a.pos.y, a.pos.z]) - x)
            if dd < best[0]: best = (dd, a.name, n)
            if dd <= 3.5 and a.element.name in ('N', 'O') and n[0] in 'NO': hb.append(f'{a.name}-{n} {dd:.2f}')
    ents = g12[(g12.pos_start <= pos) & (g12.pos_end >= pos)]
    mm = mic[mic.codon == pos]
    rows.append(dict(who_pos=pos, ecoli_pos=pos + 81, aa=rr.name, min_dist=round(best[0], 2), residue_atom=best[1], rif_atom=best[2],
        rif_moiety=moiety(best[2]), polar_contacts_le3_5A='; '.join(hb), n_grade12_entries=len(ents),
        grade12_entries=', '.join(ents.variant.str.replace('rpoB_p.', '', regex=False)), cryptic_solo_isolates=int(mm.n_solo.sum()),
        mechanism=('direct contact' if best[0] <= 4.5 else 'pocket shell' if best[0] <= 8 else 'distal/allosteric (RRDR loop)')))
at = pd.DataFrame(rows)
at.to_csv('results/rpoB_structural_atlas.tsv', sep='\t', index=False)
print(at[['who_pos', 'ecoli_pos', 'aa', 'min_dist', 'rif_moiety', 'n_grade12_entries', 'cryptic_solo_isolates', 'mechanism']].to_string(index=False))
b = mic.dropna(subset=['score']); b = b[b.n_solo >= 3]
r = spearmanr(b.score, np.log2(b.median_mic))
open('results/exploratory_distance_mic.txt', 'w').write(f'EXPLORATORY (not a locked gate). Variants with >=3 solo isolates: n={len(b)}. Spearman rho(distance, log2 median MIC) = {r.statistic:.3f}, p = {r.pvalue:.2e}\n' + b.groupby('class').agg(n=('MUTATION', 'size'), median_mic=('median_mic', 'median'), mean_frac_R=('frac_R', 'mean')).to_string() + '\n')
print(open('results/exploratory_distance_mic.txt').read())
