"""Derive rpoB/RIF subsets from raw sources (see DATA_MANIFEST.md). Run from rpoB/data."""
import pandas as pd
w = pd.read_csv('who2023_master.txt', sep='\t', low_memory=False)
w[(w.gene == 'rpoB') & (w.drug.str.contains('Rif', case=False))].to_csv('who2023_rpoB_RIF.tsv', sep='\t', index=False)
e = pd.read_parquet('z_EFFECTS.parquet').reset_index()
e[(e.GENE == 'rpoB') & (e.DRUG == 'RIF')].to_csv('cryptic_rpoB_RIF_effects.tsv.gz', sep='\t', index=False)
p = pd.read_parquet('z_UKMYC_PHENOTYPES.parquet').reset_index()
p[p.DRUG == 'RIF'].to_csv('cryptic_RIF_phenotypes.tsv.gz', sep='\t', index=False)
