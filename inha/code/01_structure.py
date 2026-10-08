#!/usr/bin/env python3
"""01_structure.py - 1ZID InhA numbering map, distances to ZID adduct, RSA, G0 positive control.
Locked rule inputs only: F1 (min heavy-atom dist to ZID), F2 anchors {S94,Y158,K165}, F3 RSA.
Also emits supporting 2IDZ NADH distances for atlas annotation (not part of the gate)."""
import gzip, os, shutil, pandas as pd, numpy as np
from Bio.PDB import MMCIFParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

D = os.path.join(os.path.dirname(__file__), '..', 'data')
R = os.path.join(os.path.dirname(__file__), '..', 'results')
os.makedirs(R, exist_ok=True)

ref = ''.join(l.strip() for l in open(os.path.join(D, 'P9WGR1.fasta')) if not l.startswith('>'))
print('P9WGR1 length:', len(ref))

p = MMCIFParser(QUIET=True)
def ungz(name):
    out = os.path.join('/tmp', name.replace('.gz',''))
    with gzip.open(os.path.join(D, name),'rb') as f, open(out,'wb') as g: shutil.copyfileobj(f,g)
    return out
s = p.get_structure('1zid', ungz('1ZID.cif.gz'))
chain = s[0]['A']
residues = [r for r in chain if is_aa(r, standard=False)]
print('chain A resolved residues:', len(residues), 'numbering', residues[0].id[1], '-', residues[-1].id[1])

# fidelity: structure sequence vs P9WGR1 (DBREF: struct 3-269 -> UNP 3-269, offset 0)
mism = []
for r in residues:
    n = r.id[1]
    a = seq1(r.resname) if r.resname != 'MSE' else 'M'
    if 1 <= n <= len(ref) and a != ref[n-1]:
        mism.append((n, ref[n-1], a, r.resname))
print('engineered/sequence mismatches vs P9WGR1:', mism)

# ZID ligand atoms (heavy)
lig = [a for r in chain for a in r if r.id[0].startswith('H_')]
zid = [a for r in chain if r.id[0] == 'H_ZID' for a in r if a.element != 'H']
print('ZID atoms:', len(zid))
zcoord = np.array([a.coord for a in zid])

rows = []
for r in residues:
    n = r.id[1]
    coords = np.array([a.coord for a in r if a.element != 'H'])
    d = np.sqrt(((coords[:, None, :] - zcoord[None, :, :])**2).sum(-1)).min()
    rows.append({'who_pos': n, 'resname': seq1(r.resname), 'min_dist': round(float(d), 3)})

# RSA on protein-only copy of chain A
import copy
prot = copy.deepcopy(chain)
for r in [r for r in prot if r.id[0] != ' ']:
    prot.detach_child(r.id)
sr = ShrakeRupley()
sr.compute(prot, level='R')
for row in rows:
    row['rsa'] = round(float(prot[row['who_pos']].sasa), 4)

df = pd.DataFrame(rows)
df['engineered'] = df.who_pos.isin([m[0] for m in mism])
df.to_csv(os.path.join(R, 'residue_distances.tsv'), sep='\t', index=False)

# G0 positive control: S94, Y158, K165 within 6.0 A of ZID
print('\nG0 POSITIVE CONTROL (must be < 6.0 A):')
ok = True
for n in (94, 158, 165):
    dd = float(df[df.who_pos == n].min_dist.iloc[0]); ok &= dd < 6.0
    print(f'  {ref[n-1]}{n}: {dd:.2f} A')
print('G0:', 'PASS' if ok else 'FAIL')

# supporting: 2IDZ is a higher-resolution (2.0 A) wild-type InhA + ZID complex; use its ZID
# distances as a resolution cross-check for the atlas (not part of the locked gate).
s2 = p.get_structure('2idz', ungz('2IDZ.cif.gz'))
c2 = s2[0]['A']
nd = [a for r in c2 if r.id[0] == 'H_ZID' for a in r if a.element != 'H']
print('2IDZ ZID atoms:', len(nd))
if nd:
    ncoord = np.array([a.coord for a in nd])
    nmap = {}
    for r in c2:
        if is_aa(r, standard=False):
            coords = np.array([a.coord for a in r if a.element != 'H'])
            nmap[r.id[1]] = round(float(np.sqrt(((coords[:, None, :] - ncoord[None, :, :])**2).sum(-1)).min()), 3)
    df['dist_zid_2idz'] = df.who_pos.map(nmap)
    df.to_csv(os.path.join(R, 'residue_distances.tsv'), sep='\t', index=False)
    same = df.dropna(subset=['dist_zid_2idz'])
    print('2IDZ cross-check: n=%d, mean |d1ZID-d2IDZ| = %.2f A, corr = %.3f' % (
        len(same), (same.min_dist - same.dist_zid_2idz).abs().mean(), same.min_dist.corr(same.dist_zid_2idz)))
