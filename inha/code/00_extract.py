#!/usr/bin/env python3
"""00_extract.py - derive inhA/fabG1 subsets from byte-locked raw sources.
Mirrors rpob/code/00_extract.py conventions. Raw inputs live in ../data and are
sha256-locked; outputs are the derived TSVs committed to the repo."""
import pandas as pd, pyarrow.parquet as pq, gzip, os

D = os.path.join(os.path.dirname(__file__), '..', 'data')

# 1. WHO 2023 master -> inhA/fabG1 x {Isoniazid, Ethionamide}
who = pd.read_csv(os.path.join(D, 'who2023_master.txt'), sep='\t', dtype=str)
sub = who[who['gene'].isin(['inhA', 'fabG1']) & who['drug'].isin(['Isoniazid', 'Ethionamide'])]
sub.to_csv(os.path.join(D, 'who2023_inhA_fabG1_INH_ETH.tsv'), sep='\t', index=False)
print('who subset rows:', len(sub))
print(sub.groupby(['drug', 'gene', 'final grade' if 'final grade' in sub.columns else 'effect']).size() if True else '')

# 2. CRyPTIC EFFECTS -> inhA/fabG1 x {INH, ETH}
eff = pq.read_table(os.path.join(D, 'z_EFFECTS.parquet')).to_pandas().reset_index()
esub = eff[eff['GENE'].isin(['inhA', 'fabG1']) & eff['DRUG'].isin(['INH', 'ETH'])]
esub.to_csv(os.path.join(D, 'cryptic_inhA_fabG1_INH_ETH_effects.tsv.gz'), sep='\t', index=False, compression='gzip')
print('cryptic effects rows:', len(esub), 'samples:', esub['UNIQUEID'].nunique())

# 3. UKMYC phenotypes -> {INH, ETH}
ph = pq.read_table(os.path.join(D, 'z_UKMYC_PHENOTYPES.parquet')).to_pandas().reset_index()
psub = ph[ph['DRUG'].isin(['INH', 'ETH'])]
psub.to_csv(os.path.join(D, 'cryptic_INH_ETH_phenotypes.tsv.gz'), sep='\t', index=False, compression='gzip')
print('phenotype rows:', len(psub))
