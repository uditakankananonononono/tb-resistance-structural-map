# rpob - rifampicin / rpoB structural atlas (builder 8)
Map WHO 2023 rpoB RIF entries and CRyPTIC clinical variants onto the Mtb RNAP-rifampicin structure (PDB 5UHC) with a pre-registered distance rule.

## Reproduce (python3, ~1 min, 2 GB RAM)
    pip install biopython pandas numpy matplotlib gemmi pyarrow scipy reportlab
    cd rpob
    # fetch the 3 gitignored raw files listed in data/DATA_MANIFEST.md into data/ and check sha256
    (cd data && python3 ../code/00_extract.py)   # optional: rebuilds the committed subsets
    python3 code/01_structure.py   # residue distances + numbering map, G0
    python3 code/02_validate.py    # G1 / G1b
    python3 code/03_candidates.py  # G2 (CRyPTIC solo-isolate MICs)
    python3 code/04_atlas.py       # atlas + exploratory distance-MIC
    python3 code/05_figures.py
    git log --format='%h %ad %s' --date=iso -- . > /tmp/gitlog.txt && python3 code/06_paper.py

## Gates and results
- GATES_LOCKED.md (commit 4cce4ec, before data), GATES_ADDENDUM_1.md (2e8dd5f, before G1)
- G0 PASS; G1 PASS 123/136 (90.4%, CI 84.3-94.3%, thin margin); G2 no candidates (4 thin, n=1)
- Paper: paper/rpoB_structural_atlas_paper.pdf
- Payload checksums: MANIFEST.sha256 (this directory)
