# katG-STRUCT-ATLAS — LOCKED SUCCESS GATES (builder 6 slice)
Locked: 2026-09-23 13:29 IST, BEFORE any WHO-catalogue agreement is measured.
Project: TB drug resistance — tb-resistance-structural-map. Slice: katG / isoniazid.

## Question
Can mutation-impact labels derived ONLY from structural analysis of PDB 1SJ2
(M. tuberculosis catalase-peroxidase KatG, 2.41 A) reproduce the WHO 2023
mutation-catalogue isoniazid-resistance grades for katG, and what is the
per-mutation structural mechanism atlas?

## Prior-art verdict: CROWDED (not DONE)
Closest works:
1. Torres et al. "Modeling the Structural Origins of Drug Resistance to Isoniazid
   via key Mutations in KatG" https://pmc.ncbi.nlm.nih.gov/articles/PMC7330162/ —
   small hand-picked mutation set, docking-based, no catalogue-scale validation.
2. Rabha et al. 2019, Microb Pathog 129:152-160
   https://www.sciencedirect.com/science/article/abs/pii/S0882401018314955 —
   double mutants only, docking+MD.
3. Saravanan et al. 2019, Sci Rep 9:6928
   https://www.nature.com/articles/s41598-019-46756-x — 98 Chennai isolates,
   SDM stability scores, no WHO-grade validation gate.
Distinct angle here: catalogue-scale (every graded katG entry), validation-FIRST
(locked >=90% concordance gate vs WHO grade-1/2 before any novel interpretation),
label-free structural scoring, then novel candidates from clinical isolates.

## Data (byte-lock in milestone 2)
- WHO 2023 catalogue v2 master file: GTB-tbsequencing/mutation-catalogue-2023
  GitHub "Final Result Files/WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
- Reference KatG: UniProt P9WIE5 (Rv1908c, H37Rv), 740 aa (verify).
- Structure: PDB 1SJ2 (2.41 A, native Mtb KatG dimer).
- Novel candidates: BV-BRC / CRyPTIC public katG variant tables.

## Structure-sequence fidelity (G2)
- Align 1SJ2 SEQRES/ATOM sequence to H37Rv katG; produce an explicit
  PDB-residue <-> gene-residue numbering map. Any engineered construct
  mutations documented; those sites excluded from atlas claims.

## Label rules (LOCKED, mechanism-based, no outcome peeking)
Computed per protein-altering katG variant:
Features: (F1) min side-chain heavy-atom distance to heme; (F2) min distance to
catalytic/adduct residues {R104,W107,H108,M255,Y229,H270,W321,D381} (numbering
verified in G2); (F3) relative solvent accessibility (Shrake-Rupley, in-house);
(F4) dimer-interface proximity (min inter-chain heavy-atom distance);
(F5) swap class: charge gain/loss, Pro gain/loss, volume change, BLOSUM62.

STRUCT-R (predicted resistance-associated) if ANY of:
- R0: loss-of-function (nonsense, frameshift, start-loss, in-frame del/ins
  touching R1-R4 sites) — mirrors WHO additional grading rule for katG LoF.
- R1: mutated residue IS a catalytic/adduct/interface-heme residue (F2 list).
- R2: any side-chain heavy atom of the WT residue within 5.0 A of heme (F1<=5.0).
- R3: buried (RSA<0.20) AND destabilizing swap (net charge change, Pro gain/loss,
  or |volume change|>80 A^3).
- R4: dimer interface (F4<=5.0 A) AND non-conservative swap (BLOSUM62<0).
Else STRUCT-neutral.

## Gates
- G1 DATA: catalogue parsed; every katG protein-altering entry extracted with
  grade; counts reported (total, per grade). Checksums logged.
- G2 FIDELITY: numbering map verified; construct anomalies documented.
- G3 POSITIVE CONTROLS (pipeline recovers known answers): S315T labeled
  STRUCT-R; >=4 of {W107,H108,R104,Y229,M255,W321,H270,D381}-site catalogue
  mutations labeled STRUCT-R; >=2 grade-4/5 neutral polymorphisms labeled
  STRUCT-neutral.
- G4 PRIMARY VALIDATION GATE: on katG missense entries with WHO grade 1/2,
  STRUCT label agreement >=90%. Grade 4/5 entries = negative set (specificity
  reported honestly, not gated). Grade 3 = uncertain: excluded from the
  denominator, reported separately. Every discrepancy documented individually.
- G5 NOVELTY (only after G4 passes): clinical-isolate katG variants absent from
  the WHO catalogue, scored by the SAME locked rules; reported as candidates.
- If G4 fails: report the failure with the full discrepancy table; any
  refinement is labeled POST-HOC EXPLORATORY, never merged into the gate.

## Discipline
Real data, real numbers. verified/thin/missing counts in every report.
Negative results preserved. No re-fishing.
