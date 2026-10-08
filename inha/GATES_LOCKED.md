# inhA slice (builder 7) - locked gates
Locked: 2026-10-08 19:01 IST, before any outcome data (WHO catalogue inhA/fabG1 grades, CRyPTIC INH/ETH phenotypes) was downloaded or inspected by this builder.

## CONTAMINATION DISCLOSURE (mandatory, on the record)
The lane record carries a contamination note: the B07 v2 lane wrote on 2026-09-23 ~6:43 PM IST
"even though a previous builder already peeked" - an earlier B07 builder may have seen WHO inhA
outcome grades before any gate was locked. This builder (v-current) states for the record:
1. No WHO 2023 catalogue row for inhA or fabG1, and no CRyPTIC INH/ETH phenotype, has been
   downloaded, opened, grepped, or otherwise inspected in this lane before the commit of this file.
2. The only facts gathered pre-lock are method inputs: PDB header metadata for 1ZID/2IDZ/4TZK
   (titles, resolution, chain range, ligand identity ZID = isonicotinic-acetyl-NAD adduct,
   DBREF to UniProt P0A5Y6), the UniProt accession lookup, and prior-art URLs. None of these
   contain WHO grades or phenotype outcomes.
3. Proof of lock-before-results is the git commit order: this file's commit precedes the data
   byte-lock commit, which precedes every results commit. Commit SHAs are reported per milestone
   and the top-level manifest records payload hashes at seal.

## Question
Can a fixed, pre-specified structural rule on the M. tuberculosis InhA - isonicotinic-acyl-NAD
(INH-NAD adduct) complex (PDB 1ZID) reproduce the WHO 2023 grade-1/2 resistance assignments for
inhA coding variants (isoniazid and ethionamide), and which grade-3 or uncatalogued clinical
inhA variants does the same rule flag as candidate resistance mutations? The fabG1 (inhA
promoter) upstream variants are handled by a separate pre-specified regulatory rule, because no
protein structure can score a promoter nucleotide.

## Prior art verdict: CROWDED (niche remains)
- PLOS One 2015, MD mechanistic study of INH resistance via InhA mutants:
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0144635 - hand-picked mutants, MD.
- Biophysical Journal 2005, MD of I21V/I16T InhA-NADH:
  https://www.cell.com/biophysj/fulltext/S0006-3495(05)72739-9 - two mutants, affinity focus.
- J Mol Graph Model 2020, integrated computational InhA mutation study:
  https://www.sciencedirect.com/science/article/abs/pii/S109332632030557X - curated set, no catalogue gate.
- Thermodynamic integration + MD INH/ETH cross-resistance (prodrug activation focus):
  https://research.mpu.edu.mo/en/publications/thermodynamic-integration-combined-with-molecular-dynamic-simulat/
- In silico functional implications of DR mutations in Mtb: https://pmc.ncbi.nlm.nih.gov/articles/PMC12703862/
- WHO mutation catalogue: 2021 1st ed. (Lancet Microbe, PIIS2666-5247(21)00301-3) and the 2023
  2nd ed. (WHO-UCN-TB-2023.6) used as the answer key by the katg/rpob slices of this repo.
Our niche: catalogue-scale (every graded inhA coding entry for INH and ETH), validation-FIRST with
a locked >=90% concordance gate before any novel interpretation, a transparent distance/contact +
stability rule with no fitted parameters, then novel candidates from CRyPTIC isolates with MIC
evidence, plus a per-mutation mechanism atlas.

## Data (sources, to be byte-locked at milestone 2, AFTER this gate commit)
- Structure (primary): PDB 1ZID, Mtb InhA + INH-NAD adduct (ligand ZID, chain A, residues 3-269,
  2.70 A; engineered construct - mutation(s) to be documented at G2-fidelity and those sites
  excluded from atlas claims).
- Structure (supporting): PDB 2IDZ (wild-type InhA + NADH, 2.0 A) for the cofactor reference;
  PDB 4TZK (1.62 A inhibitor complex) for cross-checks only.
- Reference sequence: UniProt P0A5Y6 (INHA_MYCTU, Rv1484, 269 aa), the DBREF target of 1ZID.
- Answer key: WHO 2023 catalogue, 2nd ed. master file (WHO-UCN-TB-2023.6-eng), rows with gene in
  {inhA, fabG1} and drug in {Isoniazid, Ethionamide} - same master file byte-locked by the rpob slice.
- Novel candidates: CRyPTIC public tables (Zenodo record 15680920: EFFECTS.parquet,
  UKMYC_PHENOTYPES.parquet) filtered to inhA/fabG1 with INH/ETH MICs; EBI GENOMES.csv for lineage.

## Structure-sequence fidelity (G2-fidelity)
- Align 1ZID chain A SEQRES/ATOM sequence to P0A5Y6; produce an explicit PDB-residue <->
  gene-residue numbering map. Document the engineered construct mutation(s); exclude those sites
  from atlas claims. H37Rv/WHO numbering is the reporting numbering.

## Locked rule (structural classifier, no fitted parameters)
Computed per inhA coding protein-altering variant (missense, in-frame indel, nonsense):
Features: (F1) min heavy-atom distance of the WT residue to the ZID adduct in 1ZID;
(F2) catalytic/adduct-anchor membership: {S94, Y158, K165} (fixed from literature: catalytic Tyr
and Lys, Ser94 adduct anchor); (F3) relative solvent accessibility (Shrake-Rupley, in-house, on
the 1ZID monomer); (F4) swap class: net charge change, Pro gain/loss, volume change, BLOSUM62.
STRUCT-R (predicted resistance-associated) if ANY of:
- R1: F1 <= 8.0 A (direct adduct-proximity rule; thresholds: direct contact <= 4.5 A, pocket
  shell 4.5-8.0 A - mechanism classes reported per mutation).
- R2: residue in the catalytic/adduct-anchor set {S94, Y158, K165}.
- R3: buried (RSA < 0.20) AND destabilizing swap (net charge change, Pro gain/loss, or
  |volume change| > 80 A^3).
- R0 (LoF mirror): nonsense, frameshift, start-loss, or in-frame del/ins touching R1/R2 sites -
  mirrors the WHO additional grading logic for target LoF; counted as STRUCT-R.
Else STRUCT-neutral.
Multi-residue in-frame events score min distance over the span; insertions score the two flanking
residues (mirrors the rpob addendum convention).
Promoter rule (R-PRO, separate): variants upstream of the fabG1 start codon (the fabG1-inhA
operon promoter, e.g., c.-15C>T) cannot be scored structurally. They are pre-classified as
REGULATORY-R candidates (overexpression mechanism) and are EXCLUDED from the G1 structural
denominator; their WHO grades are reported separately as a descriptive concordance line, never
merged into G1. No outcome data was used to choose this rule - the promoter position class
(-1..-50 relative to the fabG1 start) is fixed here.

## Gates
- G0 POSITIVE CONTROL (pipeline correctness; computed at analysis time, before G1): after the
  numbering map, S94, Y158 and K165 must each lie within 6.0 A of the ZID adduct. If not, the
  mapping is wrong; fix mapping only, never thresholds.
- G1 PRIMARY VALIDATION GATE: >= 90% of unique WHO 2023 grade-1/2 inhA coding protein-altering
  entries (INH and ETH, union of per-drug grade-1/2 variant sets; per-drug concordance also
  reported) must be labeled STRUCT-R. Denominator = all such entries; entries at residues
  unresolved in 1ZID count as failures, not dropped. Grade 3 = uncertain: excluded from the
  denominator, reported separately. Novel interpretation (G3) proceeds only if G1 passes.
- G1b SPECIFICITY CHECK (reported, not pass/fail): fraction of WHO grade-4/5 inhA coding entries
  labeled STRUCT-R; the rule is uninformative if this is not clearly lower than the G1 rate.
- G2 NOVELTY (only after G1 passes): grade-3 or uncatalogued inhA coding variants labeled
  STRUCT-R AND with median INH or ETH MIC above the CRyPTIC ECOFF in isolates carrying the
  variant as the sole inhA/fabG1 non-synonymous/promoter variant, n >= 3 isolates. Candidates
  failing n >= 3 are reported as thin, not as findings.
- DISCREPANCIES: every grade-1/2 entry labeled STRUCT-neutral is listed individually with a
  documented explanation attempt; none are removed.
- IF G1 FAILS: report the failure with the full discrepancy table. Any refinement is labeled
  POST-HOC EXPLORATORY, never merged into the gate. Thresholds are not re-fished.

## Discipline
Real data, real numbers; verified/thin/missing counts in every report. Negative results
preserved. Public data only (PDB, UniProt, GTB-tbsequencing GitHub, Zenodo, EBI): stop and
report if anything needs a login, a key beyond git, or payment.
