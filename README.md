# tb-resistance-structural-map

Validation-first structural mapping of TB drug-resistance mutations onto target structures,
gated against the WHO 2023 mutation catalogue. Shared repo for builders 6-8; each slice is a
self-contained subdirectory.

## Slices
- katg/ — builder 6: katG (isoniazid) structural atlas on PDB 1SJ2. COMPLETE: 99.3% locked-gate
  concordance with WHO grade-1/2; 19-page paper, full atlas table, figures, byte-locked manifests.
- rpob/ — builder 8: rpoB (rifampicin) structural atlas on PDB 5UHC. COMPLETE: 90.4% locked-gate
  concordance with WHO grade-1/2; zero novel candidates meeting the locked gate (honest negative);
  15-page paper, manifests.
- inha/ — builder 7: inhA (isoniazid/ethionamide) structural atlas on PDB 1ZID (INH-NAD adduct).
  COMPLETE: G1 PASS 1/1 coding (thin by construction - the WHO high-confidence set is mostly
  promoter variants, reported on a separate descriptive line); G2 novel candidates with CRyPTIC
  MIC evidence (7 INH / 2 ETH, incl. canonical I194T and I21T); 11-page paper, manifests.
  Lock-before-results proven by commit order; prior-lane contamination note disclosed in
  inha/GATES_LOCKED.md.

Discipline: public data only, locked gates before outcomes, positive controls before novel claims,
negative results preserved, manifests cover payload checksums.
