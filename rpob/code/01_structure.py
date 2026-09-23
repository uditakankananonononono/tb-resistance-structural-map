"""Per-residue min heavy-atom distance from RpoB (5UHC chain bound to RFP) to rifampicin.
Maps structure numbering to WHO numbering (H37Rv RpoB, 1172-aa annotation) by sequence alignment."""
import gemmi, numpy as np, pandas as pd
from Bio import Align
st = gemmi.read_structure('data/5UHC.cif.gz'); st.remove_hydrogens(); m = st[0]
rif = [a for ch in m for r in ch if r.name == 'RFP' for a in r]
rifxyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in rif])
rifchain = [ch.name for ch in m for r in ch if r.name == 'RFP'][0]
# chain containing the RpoB polymer: the chain with most residues within 5 A of RFP
best = None
for ch in m:
    pol = [r for r in ch if gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
    if len(pol) < 500: continue
    rows = []
    for r in pol:
        xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r])
        d = np.sqrt(((xyz[:, None, :] - rifxyz[None]) ** 2).sum(-1)).min()
        rows.append((r.seqid.num, r.name, gemmi.find_tabulated_residue(r.name).one_letter_code.upper(), d))
    n5 = sum(1 for x in rows if x[3] <= 5)
    if best is None or n5 > best[0]: best = (n5, ch.name, rows)
n5, chain, rows = best
df = pd.DataFrame(rows, columns=['auth_num', 'resn', 'aa', 'min_dist'])
# reference: UniProt P9WGY9 (1178 aa). WHO catalogue uses the 1172-aa annotation = UniProt minus first 6 residues.
ref = ''.join(l.strip() for l in open('data/P9WGY9.fasta') if not l.startswith('>'))
who_seq = ref[6:]
al = Align.PairwiseAligner(mode='global', open_gap_score=-10, extend_gap_score=-0.5, match_score=2, mismatch_score=-1)
struct_seq = ''.join(df.aa)
a = al.align(struct_seq, who_seq)[0]
mp = {}
for (s0, s1), (t0, t1) in zip(*a.aligned):
    for k in range(s1 - s0): mp[s0 + k] = t0 + k + 1
df['who_pos'] = [mp.get(i) for i in range(len(df))]
df['who_aa'] = [who_seq[p - 1] if p else None for p in df.who_pos]
df['chain'] = chain
df.to_csv('results/residue_distances.tsv', sep='\t', index=False)
print('RFP chain', rifchain, 'protein chain', chain, 'residues', len(df), 'n<=5A', n5)
print('identity of mapped', (df.aa == df.who_aa).mean())
print(df[df.min_dist <= 4.5][['auth_num', 'aa', 'who_pos', 'who_aa', 'min_dist']].to_string(index=False))
