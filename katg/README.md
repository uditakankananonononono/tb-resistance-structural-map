# katG structural atlas (builder 6 slice)

Validation-first mapping of clinical katG mutations onto the KatG crystal structure (PDB 1SJ2),
gated against the WHO 2023 mutation catalogue before any novel interpretation.

## Result headline
Locked mechanism-based structural rules reproduce WHO grade-1/2 katG resistance labels at
134/135 = 99.3% (gate >=90%); the frozen rules then triage 328 of 832 uncertain (grade-3)
missense variants as structurally resistance-like and nominate 16 novel candidates observed
in clinical isolates but absent from the catalogue (incl. adduct-residue p.Tyr229His).

## Layout
- code/    — pipeline (parse -> features+rules -> atlas -> figures)
- data/    — manifests + checksums ONLY (bulk data fetched from public sources)
- results/ — katg_structural_atlas.tsv (1,654 variants), labeled tables, figures (PDF)
- paper/   — katg_paper.tex + katg_paper.pdf (19 pp)
- GATES_LOCKED.md — pre-registered gates (locked before outcome data was touched)

## Reproduce
1. Fetch inputs listed in data/README.md; verify against data/byte_lock_v1.txt.
2. python3 code/parse_catalogue.py
3. python3 code/features.py            # prints G3/G4 gate results
4. python3 code/atlas_and_figures.py   # atlas + figs 1-3
5. python3 code/figures_rest.py        # figs 2,4,5,6 + stats
6. cd paper && pdflatex katg_paper.tex && pdflatex katg_paper.tex

Requires python3.10+, biopython, numpy, pandas, scipy, matplotlib. No MD/docking/GPU.
