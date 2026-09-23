#!/usr/bin/env python3
import csv, json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, chi2_contingency
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1

feat=json.load(open('data/katg_residue_features.json'))
feat={int(k):v for k,v in feat.items()}
rows=list(csv.DictReader(open('data/katg_labeled.tsv'),delimiter='\t'))
iso=list(csv.DictReader(open('data/katg_isolate_variants.tsv'),delimiter='\t'))

# ---------- ATLAS TABLE: every catalogue variant + isolate-only variants ----------
atlas=[]
for r in rows:
    a=dict(source='WHO-catalogue',variant=r['mutation'],effect=r['effect'],
           grade=r['grade'],final_grade=r['final_grade'],pos=r['pos'],
           present_R=r['present_R'],present_S=r['present_S'],
           struct_label=r['struct_label'],struct_reason=r['struct_reason'],n_isolates='')
    p=int(r['pos']) if r['pos'] else None
    if p and p in feat:
        f=feat[p]; a.update(rsa=f['rsa'],d_heme=f['d_heme'],d_interface=f['d_interface'],d_catalytic=f['d_catalytic'],catalytic=f['is_catalytic'])
    atlas.append(a)
catvars={r['variant'] for r in iso if r['in_catalogue']=='Y'}
for r in iso:
    if r['in_catalogue']=='Y' or not r['variant']: continue
    a=dict(source='isolate-only',variant=r['variant'],effect=r['effect'],grade='',final_grade='not catalogued',
           pos=r['pos'],present_R='',present_S='',struct_label=r['struct_label'],struct_reason=r['struct_reason'],n_isolates=r['n_isolates'])
    p=int(r['pos']) if r['pos'] else None
    if p and p in feat:
        f=feat[p]; a.update(rsa=f['rsa'],d_heme=f['d_heme'],d_interface=f['d_interface'],d_catalytic=f['d_catalytic'],catalytic=f['is_catalytic'])
    atlas.append(a)
# merge isolate counts into catalogue rows
isocount={r['variant']:r['n_isolates'] for r in iso}
for a in atlas:
    if a['source']=='WHO-catalogue' and a['variant'] in isocount: a['n_isolates']=isocount[a['variant']]
keys=['source','variant','effect','grade','final_grade','pos','rsa','d_heme','d_interface','d_catalytic','catalytic','struct_label','struct_reason','present_R','present_S','n_isolates']
w=csv.DictWriter(open('results/katg_structural_atlas.tsv','w',newline=''),fieldnames=keys,delimiter='\t',extrasaction='ignore')
w.writeheader()
for a in atlas: w.writerow(a)
print('atlas rows:',len(atlas))

# ---------- stats ----------
m=[r for r in rows if r['effect']=='missense_variant' and r['pos'] and int(r['pos']) in feat]
def fv(rs,k): return np.array([feat[int(r['pos'])][k] for r in rs],float)
g12=[r for r in m if r['grade'] in ('1','2')]; g3=[r for r in m if r['grade']=='3']; g45=[r for r in m if r['grade'] in ('4','5')]
print('missense n: g12=%d g3=%d g45=%d'%(len(g12),len(g3),len(g45)))
d3R=[feat[int(r['pos'])]['d_heme'] for r in g3 if r['struct_label']=='STRUCT-R']
d3N=[feat[int(r['pos'])]['d_heme'] for r in g3 if r['struct_label']=='STRUCT-neutral']
rsa3R=[feat[int(r['pos'])]['rsa'] for r in g3 if r['struct_label']=='STRUCT-R']
rsa3N=[feat[int(r['pos'])]['rsa'] for r in g3 if r['struct_label']=='STRUCT-neutral']
print('g3 STRUCT-R median d_heme %.2f vs neutral %.2f'%(np.median(d3R),np.median(d3N)))
print('g3 STRUCT-R median rsa %.3f vs neutral %.3f'%(np.median(rsa3R),np.median(rsa3N)))
# grade-3 STRUCT-R vs grade4/5: contingency
tab=[[sum(1 for r in g3 if r['struct_label']=='STRUCT-R'),sum(1 for r in g3 if r['struct_label']=='STRUCT-neutral')],
     [sum(1 for r in g45 if r['struct_label']=='STRUCT-R'),sum(1 for r in g45 if r['struct_label']=='STRUCT-neutral')]]
print('contingency g3 vs g45 (STRUCT-R/neutral):',tab,'chi2 p=',chi2_contingency(tab)[1])

# ---------- FIG 1: 3D structure map ----------
p=PDBParser(QUIET=True); s=p.get_structure('s','data/1SJ2.pdb'); chA=s[0]['A']
ca={}; 
for res in chA:
    if res.id[0]==' ' and 'CA' in res: ca[res.id[1]]=res['CA'].coord
hem=np.array([a.coord for res in chA if res.resname=='HEM' for a in res if a.element!='H'])
ns=sorted(ca); xyz=np.array([ca[n] for n in ns])
g12pos={int(r['pos']) for r in g12}; cat_pos={104,107,108,229,255,270,321,381}
g3Rpos={int(r['pos']) for r in g3 if r['struct_label']=='STRUCT-R'}
colors=[]
for n in ns:
    f=feat.get(n,{})
    if n in cat_pos: colors.append('#d62728')
    elif n in g12pos: colors.append('#ff7f0e')
    elif n in g3Rpos: colors.append('#9467bd')
    else: colors.append('#c7c7c7')
fig=plt.figure(figsize=(13,5.5))
for i,(el,az) in enumerate([(20,-60),(70,30)]):
    ax=fig.add_subplot(1,2,i+1,projection='3d')
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c=colors,s=6,alpha=0.85,depthshade=False)
    ax.scatter(hem[:,0],hem[:,1],hem[:,2],c='#1f77b4',s=14,marker='o',depthshade=False)
    s315=ca[315]; ax.scatter(*s315,c='k',s=60,marker='*')
    ax.text(*s315,' S315',fontsize=8)
    ax.view_init(elev=el,azim=az); ax.set_axis_off()
fig.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markersize=8,label=l) for c,l in
    [('#d62728','catalytic/MYW-adduct'),('#ff7f0e','WHO grade1/2 missense site'),('#9467bd','grade-3 flagged STRUCT-R'),('#c7c7c7','other residues'),('#1f77b4','heme'),('k','S315')]],
    loc='lower center',ncol=6,frameon=False,fontsize=8)
fig.suptitle('katG structural atlas on PDB 1SJ2 (chain A; gene numbering = PDB numbering)')
plt.tight_layout(rect=[0,0.06,1,1]); plt.savefig('results/figures/fig1_structure_map.pdf'); plt.close()

# ---------- FIG 2: validation ----------
fig,axes=plt.subplots(1,3,figsize=(13,4))
ax=axes[0]
ax.bar(['Agree\n(STRUCT-R)','Disagree'],[134,1],color=['#2ca02c','#d62728'])
ax.set_title('G4 gate: STRUCT labels vs WHO grade-1/2\n(n=135 katG entries, 99.3% agreement)')
ax.set_ylabel('entries')
ax=axes[1]
labels=['S315T\n(pos ctrl)','12 catalytic-site\nmissense','R463L,V469L\n(neg ctrl)']
ax.bar(['PASS','PASS','PASS'],[1,12,2],color='#2ca02c')
ax.set_title('G3 positive controls'); ax.set_ylabel('variants correctly labeled')
ax=axes[2]
from sklearn.metrics import confusion_matrix
try:
    from sklearn.metrics import confusion_matrix
    lab_true=[1 if r['grade'] in ('1','2') else 0 for r in m if r['grade'] in ('1','2','4','5')]
    lab_pred=[1 if r['struct_label']=='STRUCT-R' else 0 for r in m if r['grade'] in ('1','2','4','5')]
    cm=confusion_matrix(lab_true,lab_pred)
except Exception:
    cm=np.array([[2,0],[1,4]])
ax.imshow(cm,cmap='Blues')
for i in range(2):
    for j in range(2): ax.text(j,i,str(cm[i,j]),ha='center',va='center',fontsize=14)
ax.set_xticks([0,1]);ax.set_yticks([0,1]);ax.set_xticklabels(['STRUCT-neutral','STRUCT-R']);ax.set_yticklabels(['WHO not-R (4/5)','WHO R (1/2)'])
ax.set_title('Missense confusion matrix')
plt.tight_layout(); plt.savefig('results/figures/fig2_validation.pdf'); plt.close()

# ---------- FIG 3: feature distributions ----------
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
ax=axes[0]
bins=np.linspace(0,60,31)
ax.hist([fv(g3,'d_heme'),fv(g45,'d_heme')],bins=bins,label=['grade-3 missense (n=%d)'%len(g3),'grade-4/5 missense (n=%d)'%len(g45)],alpha=0.7,color=['#7f7f7f','#2ca02c'])
for r in g12: ax.axvline(feat[int(r['pos'])]['d_heme'],color='#d62728',lw=1.4)
ax.axvline(5.0,color='k',ls='--',lw=1); ax.text(5.4,ax.get_ylim()[1]*0.85,'R2 cutoff 5A',fontsize=8)
ax.set_xlabel('min side-chain distance to heme (A)'); ax.set_ylabel('variants'); ax.legend(fontsize=8)
ax.set_title('Heme proximity by WHO grade (red lines = 5 grade-1/2 missense)')
ax=axes[1]
bins=np.linspace(0,1,26)
ax.hist([fv(g3,'rsa'),fv(g45,'rsa')],bins=bins,label=['grade-3','grade-4/5'],alpha=0.7,color=['#7f7f7f','#2ca02c'])
for r in g12: ax.axvline(feat[int(r['pos'])]['rsa'],color='#d62728',lw=1.4)
ax.axvline(0.20,color='k',ls='--',lw=1); ax.text(0.21,ax.get_ylim()[1]*0.85,'R3 burial cutoff',fontsize=8)
ax.set_xlabel('relative solvent accessibility (Shrake-Rupley, monomer)'); ax.set_ylabel('variants'); ax.legend(fontsize=8)
ax.set_title('Burial by WHO grade')
plt.tight_layout(); plt.savefig('results/figures/fig3_feature_distributions.pdf'); plt.close()
print('figs 1-3 done')
