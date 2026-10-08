#!/usr/bin/env python3
"""05_figures.py - figures 1-4."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R = os.path.join(os.path.dirname(__file__), '..', 'results')
F = os.path.join(R, 'figures'); os.makedirs(F, exist_ok=True)
at = pd.read_csv(os.path.join(R, 'inhA_structural_atlas.tsv'), sep='\t')
who = pd.read_csv(os.path.join(R, 'who_inhA_scored.tsv'), sep='\t')

# fig1: distance profile
fig, ax = plt.subplots(figsize=(10,4.2))
ax.plot(at.who_pos, at.min_dist, lw=0.9, color='#4477aa')
ax.axhline(8.0, ls='--', c='k', lw=0.8); ax.axhline(4.5, ls=':', c='k', lw=0.8)
g = who[(who.protein_level) & (who.grade<=2)]
for _, r in g.iterrows():
    ax.scatter(r.pos_start, r.score, c='red', zorder=5, s=42)
g3 = who[(who.protein_level) & (who.grade==3)]
ax.scatter(g3.pos_start, g3.score, c='olive', zorder=4, s=10, alpha=0.6)
ax.set_xlabel('InhA residue (H37Rv numbering)'); ax.set_ylabel('min dist to INH-NAD adduct (A)')
ax.set_title('InhA residue distances to the ZID adduct (1ZID); red = WHO grade-1/2 coding, olive = grade-3')
for pos, lab in [(94,'S94A'),(194,'I194T'),(21,'I21T')]:
    ax.annotate(lab, (pos, at.set_index('who_pos').min_dist.get(pos, 0)), textcoords='offset points', xytext=(6,6), fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(F,'fig1_distance_profile.png'), dpi=160); plt.close(fig)

# fig2: STRUCT-R fraction by WHO grade (coding entries)
fig, ax = plt.subplots(figsize=(5.2,3.6))
cats = ['1-2','3','4-5']; rates=[]; ns=[]
for lo, hi in [(1,2),(3,3),(4,5)]:
    s = who[(who.protein_level)&(who.grade>=lo)&(who.grade<=hi)]
    rates.append(100*(s.STRUCT_R==True).mean() if len(s) else 0); ns.append(len(s))
bars = ax.bar(cats, rates, color=['#cc6677','#ddcc77','#117733'])
for b,n in zip(bars, ns): ax.text(b.get_x()+b.get_width()/2, b.get_height()+1, f'n={n}', ha='center', fontsize=8)
ax.set_ylabel('% labeled STRUCT-R (locked rule)'); ax.set_xlabel('WHO 2023 grade (coding entries)')
ax.set_title('Locked rule selectivity across WHO grades'); ax.set_ylim(0,115)
fig.tight_layout(); fig.savefig(os.path.join(F,'fig2_rule_by_grade.png'), dpi=160); plt.close(fig)

# fig3: MIC vs distance (solo isolates, INH)
fig, axes = plt.subplots(1,2, figsize=(10,4.2), sharey=False)
for ax, drug in zip(axes, ['INH','ETH']):
    m = pd.read_csv(os.path.join(R, f'cryptic_solo_variant_mic_{drug}.tsv'), sep='\t')
    m = m[m.n_solo>=3].dropna(subset=['score'])
    c = m.status.map({'grade1':'red','grade2':'orange','grade3':'olive','grade4':'#117733','grade5':'#117733','uncatalogued':'#88ccee'})
    ax.scatter(m.score, np.log2(m.median_mic), c=c, s=26, alpha=0.85)
    ax.axvline(8.0, ls='--', c='k', lw=0.8)
    for _, r in m.iterrows():
        if r.MUTATION in ('S94A','I194T','I21T','I21V','T162S','E196A','P12R','A190S','I95V'):
            ax.annotate(r.MUTATION, (r.score, np.log2(r.median_mic)), textcoords='offset points', xytext=(4,4), fontsize=7)
    ax.set_xlabel('min dist to adduct (A)'); ax.set_title(f'{drug} (variants with >=3 solo isolates)')
axes[0].set_ylabel('log2 median MIC (mg/L)')
fig.suptitle('Drug-site distance vs resistance level in CRyPTIC solo isolates')
fig.tight_layout(); fig.savefig(os.path.join(F,'fig3_mic_vs_distance.png'), dpi=160); plt.close(fig)

# fig4: pocket map - RSA vs distance, structural neighborhoods
fig, ax = plt.subplots(figsize=(6.4,4.6))
sc = ax.scatter(at.min_dist, at.rsa, c=['#cc6677' if pd.notna(x) else '#999999' for x in at.best_grade], s=14, alpha=0.75)
ax.axvline(8.0, ls='--', c='k', lw=0.8); ax.axhline(0.20, ls=':', c='k', lw=0.8)
ax.set_xlabel('min dist to INH-NAD adduct (A)'); ax.set_ylabel('relative solvent accessibility')
ax.set_title('InhA structural neighborhoods (red = residue carries a WHO-graded coding entry)')
for pos in (94, 158, 165, 194, 21, 16, 95):
    row = at[at.who_pos==pos]
    if len(row): ax.annotate(f'{row.resname.iloc[0]}{pos}', (row.min_dist.iloc[0], row.rsa.iloc[0]), textcoords='offset points', xytext=(5,5), fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(F,'fig4_pocket_map.png'), dpi=160); plt.close(fig)
print('figures written:', sorted(os.listdir(F)))
