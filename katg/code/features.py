#!/usr/bin/env python3
"""katG structural atlas: per-residue structural features from PDB 1SJ2, locked rules -> STRUCT labels."""
import csv, re, sys, numpy as np
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1

CATALYTIC={104,107,108,229,255,270,321,381}
MAXASA={'A':121,'R':265,'N':187,'D':187,'C':148,'Q':214,'E':214,'G':97,'H':216,'I':195,'L':191,'K':230,'M':203,'F':228,'P':154,'S':143,'T':163,'W':264,'Y':255,'V':165}
VOL={'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
CHARGE={'R':1,'K':1,'H':0.5,'D':-1,'E':-1}
RADIUS={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
B62={}
def load_blosum():
    global B62
    try:
        from Bio.Align import substitution_matrices
        m=substitution_matrices.load("BLOSUM62")
        aa=list(m.alphabet)
        for a in aa:
            for b in aa: B62[(a,b)]=m[a,b]
    except Exception as e:
        print('BLOSUM62 load failed:',e); sys.exit(1)
load_blosum()

p=PDBParser(QUIET=True)
s=p.get_structure('s','data/1SJ2.pdb')
chA=s[0]['A']; chB=s[0]['B']
def heavy_atoms(chain, sidechain_only=False, protein_only=True):
    out=[]
    for res in chain:
        if res.id[0]!=' ': continue
        for at in res:
            if at.element=='H': continue
            if sidechain_only and at.name in ('N','CA','C','O'): continue
            out.append((res.id[1],at))
    return out
# protein atoms per chain
atomsA=[at for res in chA if res.id[0]==' ' for at in res if at.element!='H']
atomsB=[at for res in chB if res.id[0]==' ' for at in res if at.element!='H']
coordsB=np.array([a.coord for a in atomsB])
# heme atoms (chain A HEM 1500)
hem=[at for res in chA if res.resname=='HEM' for at in res if at.element!='H']
hemC=np.array([a.coord for a in hem])
print('heme atoms:',len(hem),'| chainA atoms:',len(atomsA),'| chainB atoms:',len(atomsB))

# Shrake-Rupley SASA on chain A (monomer context) 
def shrake_rupley(atoms, probe=1.4, npts=100):
    coords=np.array([a.coord for a in atoms])
    radii=np.array([RADIUS.get(a.element,1.70)+probe for a in atoms])
    # neighbor cutoff via grid-free broadcasting in chunks
    sasas=np.zeros(len(atoms))
    # golden spiral sphere
    idx=np.arange(npts); phi=np.pi*(3-np.sqrt(5))
    y=1-2*(idx+0.5)/npts; r=np.sqrt(1-y*y); t=phi*idx
    sph=np.stack([np.cos(t)*r,y,np.sin(t)*r],axis=1)
    N=len(atoms)
    for i in range(N):
        d=np.linalg.norm(coords-coords[i],axis=1)
        nb=np.where((d<radii+radii[i])&(d>0))[0]
        pts=coords[i]+radii[i]*sph
        acc=np.ones(npts,bool)
        for j in nb:
            if not acc.any(): break
            dd=np.linalg.norm(pts-coords[j],axis=1)
            acc &= dd>=radii[j]-1e-6
        sasas[i]=4*np.pi*radii[i]**2*acc.sum()/npts
    return sasas
res_list=[res for res in chA if res.id[0]==' ']
atom_owner=[]
for i,a in enumerate(atomsA): atom_owner.append(a.get_parent().id[1])
sasa=shrake_rupley(atomsA)
res_sasa={}
for i,a in enumerate(atomsA):
    rn=atom_owner[i]; res_sasa[rn]=res_sasa.get(rn,0)+sasa[i]

def res_atoms(resnum,chain=chA):
    try: res=chain[(' ',resnum,' ')]
    except KeyError: return []
    return [at for at in res if at.element!='H' and at.name not in ('N','CA','C','O')]
feat={}
for res in res_list:
    n=res.id[1]; rn=res.resname
    aa1=protein_letters_3to1.get(rn.capitalize(),'?')
    sc=[at for at in res if at.element!='H' and at.name not in ('N','CA','C','O')]
    if not sc:  # glycine: use CA
        sc=[at for at in res if at.name=='CA']
    scC=np.array([a.coord for a in sc])
    d_hem=float(np.min(np.linalg.norm(hemC[None,:,:]-scC[:,None,:],axis=2)))
    d_if=float(np.min(np.linalg.norm(coordsB[None,:,:]-scC[:,None,:],axis=2)))
    d_cat=min((float(np.min(np.linalg.norm(np.array([a.coord for a in res_atoms(c)])[None,:,:]-scC[:,None,:],axis=2))) for c in CATALYTIC if c!=n and res_atoms(c)), default=np.nan)
    rsa=res_sasa.get(n,0)/MAXASA.get(aa1,200)
    feat[n]=dict(pos=n,wt=aa1,rsa=round(rsa,4),d_heme=round(d_hem,2),d_interface=round(d_if,2),d_catalytic=round(d_cat,2) if d_cat==d_cat else '',is_catalytic=(n in CATALYTIC))

print('residues with features:',len(feat))
import json
json.dump(feat,open('data/katg_residue_features.json','w'))

# apply locked rules to catalogue
def struct_label(wt,pos,mt,effect):
    if effect in ('frameshift','stop_gained','start_lost','LoF','feature_ablation','initiator_codon_variant'): return 'STRUCT-R','R0-LoF'
    if effect!='missense_variant' or not pos: return 'NA','non-missense'
    f=feat.get(int(pos))
    if f is None: return 'NA','no-structure(res 2-23 disordered)'
    reasons=[]
    if f['is_catalytic']: reasons.append('R1-catalytic/adduct site')
    if f['d_heme']<=5.0: reasons.append('R2-heme<=5A')
    cwt,cm=CHARGE.get(wt,0),CHARGE.get(mt,0)
    chg=(np.sign(cwt)!=np.sign(cm)) or (cwt==0 and cm!=0) or (cwt!=0 and cm==0)
    pro=(wt=='P')!=(mt=='P')
    dv=abs(VOL.get(mt,0)-VOL.get(wt,0))
    if f['rsa']<0.20 and (chg or pro or dv>80): reasons.append(f'R3-buried+destab(chg={bool(chg)},pro={pro},dV={dv:.0f})')
    b=B62.get((wt,mt),B62.get((mt,wt),0))
    if f['d_interface']<=5.0 and b<0: reasons.append('R4-interface+nonconservative')
    return ('STRUCT-R' if reasons else 'STRUCT-neutral'),';'.join(reasons) if reasons else 'none'

rows=list(csv.DictReader(open('data/katg_catalogue.tsv'),delimiter='\t'))
labeled=0
for r in rows:
    lab,why=struct_label(r['wt_aa'],r['pos'],r['mut_aa'],r['effect'])
    r['struct_label']=lab; r['struct_reason']=why
    if lab!='NA': labeled+=1
w=csv.DictWriter(open('data/katg_labeled.tsv','w',newline=''),fieldnames=list(rows[0].keys()),delimiter='\t')
w.writeheader(); w.writerows(rows)
print('labeled (structurally addressable):',labeled,'/',len(rows))

# G3 positive controls
by_mut={r['mutation']:r for r in rows}
print('\n--- G3 positive controls ---')
print('S315T:',by_mut['p.Ser315Thr']['struct_label'],by_mut['p.Ser315Thr']['struct_reason'])
cat_sites=['Trp107','His108','Arg104','Tyr229','Met255','Trp321','His270','Asp381']
cat_rows=[r for r in rows if r['effect']=='missense_variant' and r['pos'] and any(r['mutation'].startswith('p.'+cs) for cs in cat_sites)]
ncat_R=sum(1 for r in cat_rows if r['struct_label']=='STRUCT-R')
print('catalogue missense at catalytic sites:',len(cat_rows),'STRUCT-R:',ncat_R)
neut=[r for r in rows if r['mutation'] in ('p.Arg463Leu','p.Val469Leu')]
for r in neut: print('negative control',r['mutation'],r['struct_label'],r['struct_reason'])
