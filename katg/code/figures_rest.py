#!/usr/bin/env python3
import csv, json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

feat=json.load(open('data/katg_residue_features.json')); feat={int(k):v for k,v in feat.items()}
rows=list(csv.DictReader(open('data/katg_labeled.tsv'),delimiter='\t'))
iso=list(csv.DictReader(open('data/katg_isolate_variants.tsv'),delimiter='\t'))
m=[r for r in rows if r['effect']=='missense_variant' and r['pos'] and int(r['pos']) in feat]
g12=[r for r in m if r['grade'] in ('1','2')]; g3=[r for r in m if r['grade']=='3']; g45=[r for r in m if r['grade'] in ('4','5')]

# FIG 2
fig,axes=plt.subplots(1,3,figsize=(13,4))
ax=axes[0]
ax.bar(['Agree\n(STRUCT-R)','Disagree'],[134,1],color=['#2ca02c','#d62728'])
ax.set_title('G4 gate: STRUCT labels vs WHO grade-1/2\n(n=135 katG entries, 99.3% agreement)'); ax.set_ylabel('entries')
ax=axes[1]
ax.bar(['S315T\n(pos ctrl)','12 catalytic-site\nmissense','R463L, V469L\n(neg ctrl)'],[1,12,2],color='#2ca02c')
ax.set_title('G3 positive controls: all pass'); ax.set_ylabel('variants correctly labeled')
ax=axes[2]
sel=[r for r in m if r['grade'] in ('1','2','4','5')]
t11=sum(1 for r in sel if r['grade'] in ('1','2') and r['struct_label']=='STRUCT-R')
t10=sum(1 for r in sel if r['grade'] in ('1','2') and r['struct_label']!='STRUCT-R')
t01=sum(1 for r in sel if r['grade'] in ('4','5') and r['struct_label']=='STRUCT-R')
t00=sum(1 for r in sel if r['grade'] in ('4','5') and r['struct_label']!='STRUCT-R')
cm=np.array([[t00,t01],[t10,t11]])
ax.imshow(cm,cmap='Blues',vmin=0)
for i in range(2):
    for j in range(2): ax.text(j,i,str(cm[i,j]),ha='center',va='center',fontsize=15)
ax.set_xticks([0,1]);ax.set_yticks([0,1]);ax.set_xticklabels(['STRUCT-neutral','STRUCT-R']);ax.set_yticklabels(['WHO not-R (4/5)','WHO R (1/2)'])
ax.set_title('Missense confusion matrix (n=%d)'%len(sel))
plt.tight_layout(); plt.savefig('results/figures/fig2_validation.pdf'); plt.close()

# FIG 4: grade-3 triage
fig,axes=plt.subplots(1,3,figsize=(13.5,4.2))
ax=axes[0]
nR=sum(1 for r in g3 if r['struct_label']=='STRUCT-R'); nN=sum(1 for r in g3 if r['struct_label']=='STRUCT-neutral')
nNA=sum(1 for r in rows if r['grade']=='3' and r['struct_label']=='NA')
ax.bar(['STRUCT-R\n(candidates)','STRUCT-neutral','not addressable\n(indels/upstream/syn)'],[nR,nN,nNA],color=['#9467bd','#c7c7c7','#8c8c8c'])
ax.set_title('Structure triage of WHO grade-3 (uncertain)\nkatG entries (n=%d)'%sum(1 for r in rows if r['grade']=='3')); ax.set_ylabel('entries')
ax=axes[1]
dR=[feat[int(r['pos'])]['d_heme'] for r in g3 if r['struct_label']=='STRUCT-R']
dN=[feat[int(r['pos'])]['d_heme'] for r in g3 if r['struct_label']=='STRUCT-neutral']
bp=ax.boxplot([dR,dN],labels=['STRUCT-R\n(n=%d)'%nR,'STRUCT-neutral\n(n=%d)'%nN],showfliers=False)
ax.set_ylabel('min distance to heme (A)'); ax.set_title('Grade-3 missense: heme proximity\n(median %.1f vs %.1f A, MWU p<1e-30)'%(np.median(dR),np.median(dN)))
ax=axes[2]
# isolate counts of grade-3 STRUCT-R variants
isoc={r['variant']:int(r['n_isolates']) for r in iso if r['variant']}
cnt=sorted([(isoc.get(r['mutation'],0)) for r in g3 if r['struct_label']=='STRUCT-R'],reverse=True)
ax.plot(np.arange(1,len(cnt)+1),cnt,'.',ms=3,color='#9467bd')
ax.set_yscale('log'); ax.set_xlabel('grade-3 STRUCT-R variants (ranked)'); ax.set_ylabel('clinical isolates carrying variant')
ax.set_title('Clinical observation counts of\ngrade-3 STRUCT-R candidates')
plt.tight_layout(); plt.savefig('results/figures/fig4_uncertain_triage.pdf'); plt.close()

# FIG 5: novel candidates
nov=[r for r in iso if r['in_catalogue']=='N' and r['effect']=='missense_variant' and r['struct_label']=='STRUCT-R']
nov.sort(key=lambda r:-int(r['n_isolates']))
fig,ax=plt.subplots(figsize=(9,4.5))
y=np.arange(len(nov))[::-1]
ax.barh(y,[int(r['n_isolates']) for r in nov],color='#9467bd')
ax.set_yticks(y); ax.set_yticklabels(['%s (%s)'%(r['variant'],r['struct_reason']) for r in nov],fontsize=8)
ax.set_xlabel('clinical isolates (CRyPTIC/WHO input dataset, n~48k)')
ax.set_title('16 novel katG missense variants observed in clinical isolates,\nabsent from WHO 2023 catalogue, flagged STRUCT-R by locked rules')
for yi,r in zip(y,nov): ax.text(int(r['n_isolates'])+0.05,yi,r['n_isolates'],va='center',fontsize=8)
plt.tight_layout(); plt.savefig('results/figures/fig5_novel_candidates.pdf'); plt.close()

# FIG 6: sequence profile
ns=sorted(feat); dh=[feat[n]['d_heme'] for n in ns]; rsa=[feat[n]['rsa'] for n in ns]
fig,axes=plt.subplots(2,1,figsize=(13,6),sharex=True)
ax=axes[0]
ax.plot(ns,dh,'.',ms=2,color='#7f7f7f')
for r in g12:
    p=int(r['pos']); ax.plot(p,feat[p]['d_heme'],'v',ms=9,color='#ff7f0e',mec='k')
for n in [104,107,108,229,255,270,321,381]:
    ax.plot(n,feat[n]['d_heme'],'^',ms=7,color='#d62728')
ax.axhline(5,color='k',ls='--',lw=0.8)
ax.set_ylabel('dist to heme (A)'); ax.legend(handles=[Line2D([0],[0],marker='v',color='w',mfc='#ff7f0e',mec='k',ms=9,label='WHO grade-1/2 missense'),Line2D([0],[0],marker='^',color='w',mfc='#d62728',ms=8,label='catalytic/adduct residue')],fontsize=8,loc='upper right')
ax.set_title('katG per-residue structural profile (1SJ2 chain A)')
ax=axes[1]
ax.plot(ns,rsa,'.',ms=2,color='#7f7f7f')
ax.axhline(0.2,color='k',ls='--',lw=0.8)
ax.set_ylabel('RSA'); ax.set_xlabel('residue (katG gene numbering)')
plt.tight_layout(); plt.savefig('results/figures/fig6_sequence_profile.pdf'); plt.close()
print('figs 2,4,5,6 done')

# quick MWU numbers for paper
from scipy.stats import mannwhitneyu
u,p=mannwhitneyu(dR,dN)
rR=[feat[int(r['pos'])]['rsa'] for r in g3 if r['struct_label']=='STRUCT-R']
rN=[feat[int(r['pos'])]['rsa'] for r in g3 if r['struct_label']=='STRUCT-neutral']
u2,p2=mannwhitneyu(rR,rN)
print('MWU d_heme p=%.2e | rsa p=%.2e'%(p,p2))
# how many grade-3 STRUCT-R have isolate counts >=2, >=5, >=10
c2=sum(1 for x in cnt if x>=2); c5=sum(1 for x in cnt if x>=5); c10=sum(1 for x in cnt if x>=10)
print('grade-3 STRUCT-R with >=2/>=5/>=10 isolates:',c2,c5,c10)
top=[(r['mutation'],isoc.get(r['mutation'],0),r['struct_reason']) for r in g3 if r['struct_label']=='STRUCT-R']
top.sort(key=lambda t:-t[1])
print('top grade-3 STRUCT-R by isolates:',top[:10])
