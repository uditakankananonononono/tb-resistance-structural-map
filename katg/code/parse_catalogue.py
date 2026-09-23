#!/usr/bin/env python3
"""Parse WHO 2023 TB mutation catalogue master file -> katG/INH table."""
import csv, re, sys, hashlib, json
from Bio.Data.IUPACData import protein_letters_3to1

IN='data/who_catalogue_master.txt'; OUT='data/katg_catalogue.tsv'
rows=[]
with open(IN) as f:
    rd=csv.reader(f,delimiter='\t')
    hdr=next(rd)
    for r in rd:
        if len(r)<106: continue
        if r[0]=='Isoniazid' and r[1]=='katG': rows.append(r)
print('katG/INH rows:',len(rows))

def parse_p(mut):
    """Return (wt,pos,mut,type) from p. notation."""
    m=re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$',mut)
    if m: return protein_letters_3to1[m.group(1)],int(m.group(2)),protein_letters_3to1[m.group(3)],'missense'
    m=re.match(r'p\.([A-Z][a-z]{2})(\d+)\*$',mut)
    if m: return protein_letters_3to1[m.group(1)],int(m.group(2)),'*','stop_gained'
    m=re.match(r'p\.([A-Z][a-z]{2})(\d+)fs$',mut)
    if m: return protein_letters_3to1[m.group(1)],int(m.group(2)),'fs','frameshift'
    m=re.match(r'p\.Met1\?$',mut)
    return None,None,None,None

out=[]
for r in rows:
    mut=r[2]; effect=r[5]; grade=r[105].strip()
    grade_num=grade[0] if grade and grade[0].isdigit() else ''
    wt,pos,mt,pt=parse_p(mut)
    if pt is None and effect in ('missense_variant','stop_gained','frameshift','start_lost','initiator_codon_variant'):
        # try looser patterns
        m=re.match(r'p\.([A-Z][a-z]{2})(\d+)',mut)
        if m: wt,pos=protein_letters_3to1.get(m.group(1),'?'),int(m.group(2)); mt='?'; pt=effect
    out.append(dict(mutation=mut,variant=r[3],tier=r[4],effect=effect,
        genomic_position=r[6],present_R=r[11],present_S=r[12],
        final_grade=grade,grade=grade_num,wt_aa=wt or '',pos=pos or '',mut_aa=mt or '',parsed_type=pt or ''))

with open(OUT,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0].keys()),delimiter='\t')
    w.writeheader(); w.writerows(out)

# summary
from collections import Counter
print('by grade:',dict(Counter(o['grade'] for o in out)))
print('by effect:',dict(Counter(o['effect'] for o in out)))
miss=[o for o in out if o['effect']=='missense_variant']
print('missense:',len(miss),'unique positions:',len(set(o['pos'] for o in miss)))
g12=[o for o in out if o['grade'] in ('1','2')]
print('grade1/2 total:',len(g12),' missense:',sum(1 for o in g12 if o['effect']=='missense_variant'),
      ' LoF-class:',sum(1 for o in g12 if o['effect'] in ('frameshift','stop_gained','start_lost','LoF','feature_ablation','initiator_codon_variant')))
g45=[o for o in out if o['grade'] in ('4','5')]
print('grade4/5 total:',len(g45),' missense:',[o['mutation'] for o in g45 if o['effect']=='missense_variant'])
unparsed=[o['mutation'] for o in out if o['effect'] in ('missense_variant','stop_gained','frameshift') and not o['pos']]
print('unparsed protein notations:',len(unparsed), unparsed[:10])
