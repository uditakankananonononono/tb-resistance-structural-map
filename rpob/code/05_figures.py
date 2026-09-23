import gemmi, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = pd.read_csv('results/residue_distances.tsv', sep='\t').dropna(subset=['who_pos'])
w = pd.read_csv('results/who_rpoB_scored.tsv', sep='\t')
mic = pd.read_csv('results/cryptic_solo_variant_mic.tsv', sep='\t')
at = pd.read_csv('results/rpoB_structural_atlas.tsv', sep='\t')
# F1 distance profile
fig, ax = plt.subplots(figsize=(10, 4))
s = d[(d.who_pos >= 150) & (d.who_pos <= 700)]
ax.plot(s.who_pos, s.min_dist, color='0.6', lw=0.8, label='min distance to RIF (5UHC)')
for g, c, mk in [(3, '#999900', '.'), (2, '#e08000', 'o'), (1, '#c00000', 'o')]:
    x = w[(w.grade == g) & w.protein_level].dropna(subset=['score'])
    x = x[(x.pos_start >= 150) & (x.pos_start <= 700)]
    ax.scatter(x.pos_start, x.score, s=14 if g < 3 else 6, c=c, marker=mk, label=f'WHO grade {g} entries', zorder=3)
ax.axhline(8, ls='--', c='k', lw=0.8); ax.axhline(4.5, ls=':', c='k', lw=0.8)
ax.axvspan(426, 452, color='#cfe8ff', alpha=0.5, label='RRDR 426-452')
ax.set_ylim(0, 60); ax.set_xlabel('RpoB residue (WHO / H37Rv numbering)'); ax.set_ylabel('Min heavy-atom distance (A)')
ax.legend(fontsize=7, loc='upper right'); ax.set_title('Figure 1. RpoB residue distance to rifampicin with WHO 2023 entries')
fig.tight_layout(); fig.savefig('results/figures/fig1_distance_profile.png', dpi=200); plt.close()
# F2 class by grade
w['gg'] = w.grade.map({1: 'Grade 1-2', 2: 'Grade 1-2', 3: 'Grade 3', 4: 'Grade 4-5', 5: 'Grade 4-5'})
t = pd.crosstab(w[w.protein_level].gg, w[w.protein_level]['class'], normalize='index')[['contact', 'shell', 'distal', 'unresolved/unscorable']]
n = w[w.protein_level].gg.value_counts()
fig, ax = plt.subplots(figsize=(6, 4))
t.plot.bar(stacked=True, ax=ax, color=['#c00000', '#f0a000', '#8090a0', '#dddddd'], rot=0)
ax.set_xticklabels([f'{i}\n(n={n[i]})' for i in t.index]); ax.set_ylabel('Fraction of entries'); ax.set_xlabel('')
ax.legend(fontsize=8, bbox_to_anchor=(1, 1)); ax.set_title('Figure 2. Structural class by WHO grade')
fig.tight_layout(); fig.savefig('results/figures/fig2_class_by_grade.png', dpi=200); plt.close()
# F3 MIC vs distance
b = mic.dropna(subset=['score']); b = b[b.n_solo >= 3]
fig, ax = plt.subplots(figsize=(7, 4.5))
cols = {'grade1': '#c00000', 'grade2': '#e08000', 'grade3': '#999900', 'grade4': '#4060c0', 'grade5': '#4060c0', 'uncatalogued': 'k'}
for st_, g in b.groupby('status'):
    ax.scatter(g.score, g.median_mic, s=12 + 4 * np.sqrt(g.n_solo), c=cols[st_], alpha=0.75, label=st_, edgecolor='w', lw=0.4)
for k, (_, r) in enumerate(b[b.n_solo >= 40].sort_values('score').iterrows()): ax.annotate(f'{r.MUTATION} (n={r.n_solo})', (r.score, r.median_mic), fontsize=6, xytext=(8 + 30 * (k % 3), 10 - 9 * (k % 4)), textcoords='offset points', arrowprops=dict(arrowstyle='-', lw=0.3))
ax.set_yscale('log', base=2); ax.axhline(0.5, ls='--', c='k', lw=0.8); ax.axvline(8, ls=':', c='k', lw=0.8)
ax.set_xlabel('Min distance of mutated residue to RIF (A)'); ax.set_ylabel('Median RIF MIC, solo isolates (mg/L)')
ax.legend(fontsize=7); ax.set_title('Figure 3. CRyPTIC MIC vs structural distance (variants with >=3 solo isolates)', fontsize=9)
fig.tight_layout(); fig.savefig('results/figures/fig3_mic_vs_distance.png', dpi=200); plt.close()
# F4 pocket projection
st = gemmi.read_structure('data/5UHC.cif.gz'); m = st[0]; ch = m[d.chain.iloc[0]]
rif = np.array([[a.pos.x, a.pos.y, a.pos.z] for c in m for r in c if r.name == 'RFP' for a in r if a.element.name != 'H'])
pk = d[d.min_dist <= 12]
ca = np.array([[r['CA'][0].pos.x, r['CA'][0].pos.y, r['CA'][0].pos.z] for r in ch if r.seqid.num in set(pk.auth_num)])
nums = [r.seqid.num - 6 for r in ch if r.seqid.num in set(pk.auth_num)]
X = np.vstack([rif, ca]); mu = X.mean(0); U, S, Vt = np.linalg.svd(X - mu); P = (X - mu) @ Vt[:2].T
cnt = dict(zip(at.who_pos, at.n_grade12_entries))
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(P[:len(rif), 0], P[:len(rif), 1], c='#2a9d2a', s=18, label='rifampicin atoms')
pc = P[len(rif):]
cc = [cnt.get(n, 0) for n in nums]
sc = ax.scatter(pc[:, 0], pc[:, 1], c=cc, cmap='Reds', s=[20 + 12 * c for c in cc], edgecolor='k', lw=0.4, vmin=0)
for (x, y), n, c in zip(pc, nums, cc):
    if c > 0 or n in (454, 493): ax.annotate(str(n), (x, y), fontsize=6, xytext=(2, 2), textcoords='offset points')
plt.colorbar(sc, label='# WHO grade-1/2 entries at residue'); ax.set_xlabel('PC1 (A)'); ax.set_ylabel('PC2 (A)'); ax.legend(fontsize=7)
ax.set_title('Figure 4. RIF pocket (residues within 12 A; CA atoms), PCA projection', fontsize=9); ax.set_aspect('equal')
fig.tight_layout(); fig.savefig('results/figures/fig4_pocket_map.png', dpi=200); plt.close()
print('ok')
