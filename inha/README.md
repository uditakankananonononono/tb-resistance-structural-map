# inhA slice (builder 7) - isoniazid/ethionamide, PDB 1ZID INH-NAD adduct

Can mutation-impact labels derived ONLY from structural analysis of InhA reproduce the WHO
2023 catalogue grade-1/2 resistance assignments for isoniazid/ethionamide?

## Results (all locked-gate, see GATES_LOCKED.md + addenda)
- G0 positive control PASS: S94 3.48 A, Y158 3.80 A, K165 2.88 A from the ZID adduct.
- G1 validation PASS, thin by construction: 1/1 unique WHO grade-1/2 coding variant (S94A,
  3.48 A direct contact). The WHO high-confidence set for this target is 14 promoter rows + 1
  coding variant - structure can only grade the coding one; promoter line reported separately.
- G1b specificity: 33.1% of grade-3 coding entries flagged, 0% of grade-4/5 coding entries.
- G2 novel candidates (locked gate: STRUCT-R, >=3 solo CRyPTIC isolates, median MIC > ECOFF):
  INH: I194T (n=34, 3.2 mg/L), T162S (n=21), E196A (n=10), P12R (n=9), I21V (n=5), I21T (n=5),
  A190S (n=3). ETH: I194T (n=34, 16 mg/L), I21T (n=5). Thin (n<3, documented): I95V, V20I,
  A93S, G40W (INH); I95V, I16T (ETH).
- Exploratory: adduct distance vs solo-isolate MIC, INH Spearman rho = -0.59, p = 3.3e-02.
- Paper: paper/inhA_structural_atlas_paper.pdf (11 pp).

## Integrity
Gates committed (8ea16ec) BEFORE data byte-lock (9591ba4) BEFORE addenda (8ea202d, bcd9baf)
BEFORE results (b9874e3). Contamination note from the lane record (a prior builder may have
seen WHO inhA grades) is disclosed in GATES_LOCKED.md; commit order is the lock-before-results
proof. The one outcome-informed correction (addendum 2, promoter-class mapping) touches only
the descriptive promoter line and is labeled as such.

## Layout
- GATES_LOCKED.md, GATES_ADDENDUM_1.md, GATES_ADDENDUM_2.md - pre-registration
- code/00_extract.py ... 06_paper.py - pipeline, in order
- data/ - byte-locked sources (raw bulk gitignored with URL+sha256; see DATA_MANIFEST.md)
- results/ - scored WHO table, atlas, CRyPTIC MIC tables, G2 outputs, figures, G1 log
- paper/ - inhA_structural_atlas_paper.pdf
