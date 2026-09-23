# rpoB slice (builder 8) - locked gates
Locked: 2026-09-23 13:40 IST, before any outcome data (WHO catalogue grades, CRyPTIC/BV-BRC phenotypes) was downloaded or inspected.

## Question
Can a fixed, pre-specified structural rule on the M. tuberculosis RNAP-rifampicin complex explain WHO 2023 grade-1/2 rpoB rifampicin-resistance mutations, and which uncertain (grade 3) or uncatalogued clinical rpoB variants does the same rule flag as candidate resistance mutations?

## Prior art verdict: CROWDED (niche remains)
- Portelli et al. 2020, Sci Rep, structure-based ML beyond RRDR: https://www.nature.com/articles/s41598-020-74648-y
- Structural basis of Mtb RpoB clinical mutations on RIF binding (Molecules 2022): https://www.mdpi.com/1420-3049/27/3/885
- Structure-informed ML for RIF resistance (bioRxiv 2024 / PMC 2025): https://doi.org/10.1101/2024.08.15.608097 , https://pmc.ncbi.nlm.nih.gov/articles/PMC12208608/
- Context: E. coli RIF binding-site deep scan, Nature 2023: https://www.nature.com/articles/s41586-023-06495-6
Our niche: a transparent, locked distance/contact rule validated against the WHO 2023 (2nd ed.) catalogue, applied to grade-3 and uncatalogued CRyPTIC variants with MIC evidence, plus a per-mutation mechanism atlas (direct contact / pocket-shaping / allosteric-distal).

## Data (sources, to be byte-locked at milestone 2)
- Structure: PDB 5UHC (Mtb RNAP initiation complex + rifampicin); fallback 5UH6 / other Mtb RNAP-RIF entries if ligand absent.
- Answer key: WHO 2023 catalogue of Mtbc mutations, 2nd edition, rpoB rows for RIF.
- Clinical isolates: CRyPTIC public tables (Zenodo) with RIF MICs; BV-BRC AMR tables if CRyPTIC unreachable.
- Sequence conservation: UniProt RpoB homologs (for tie-breaking only, not the primary rule).

## Locked rule (structural classifier)
Residue r (Mtb H37Rv numbering, P9WGY9) is "resistance-plausible" if its minimum heavy-atom distance to rifampicin in the structure is <= 8.0 A. Mechanism classes: direct contact <= 4.5 A; pocket shell 4.5-8.0 A; distal > 8.0 A. Thresholds are fixed now and will not be tuned.

## Gates
- G0 positive control (pipeline correctness): after numbering mapping, residues S450, H445 and D435 (Mtb numbering; E. coli S531, H526, D516) must each be within 6 A of rifampicin. If not, the mapping is wrong; fix mapping only, never thresholds.
- G1 validation gate: >= 90% of unique WHO 2023 grade-1/2 rpoB missense (and in-frame) entries for RIF must fall at residues classed resistance-plausible. Denominator = all such entries; entries at residues unresolved in the structure are counted as failures, not dropped. Novel interpretation proceeds only if G1 passes.
- G1b specificity check (reported, not a pass/fail gate): fraction of WHO grade-4/5 rpoB entries classed resistance-plausible; rule is uninformative if this is not lower than the G1 rate.
- G2 novel candidates: grade-3/uncatalogued rpoB variants classed resistance-plausible AND with median RIF MIC above the CRyPTIC ECOFF in isolates carrying them as the sole rpoB non-synonymous variant (n >= 3 isolates). Candidates failing n>=3 are reported as thin, not as findings.
- All discrepancies (grade-1/2 entries outside 8 A) are listed with a documented explanation attempt; none are removed.
- If G1 fails, the result is reported as negative; thresholds are not re-fished.
