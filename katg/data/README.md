# katg/data — manifests and checksums only (no bulk data)

Bulk inputs are public and fetched fresh by the pipeline (see katg/README.md):
- WHO 2023 catalogue master file (WHO-UCN-TB-2023.6-eng), GTB-tbsequencing/mutation-catalogue-2023 GitHub
- PDB 1SJ2 (files.rcsb.org/download/1SJ2.pdb)
- UniProt P9WIE5 FASTA (rest.uniprot.org)
- H37Rv NC_000962.3 katG region (NCBI eutils, complement 2153889..2156111)
- INH tier1/tier2 isolate genotype inputs (same GitHub repo, "Input data files for Solo algorithms")

byte_lock_v1.txt = SHA-256 of every input as downloaded on 2026-09-23.
