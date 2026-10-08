# Gates addendum 1 - operational definitions (written before G1 was computed)
Written 2026-10-08 ~19:08 IST. At this point the WHO subset has been extracted and its column
layout and effect-type counts viewed (grade distribution by count only); no distances have been
joined to WHO entries and G1 has not been computed. Thresholds and gates are unchanged; this
only states how entries are scored.

1. Grade parsing: WHO column "FINAL CONFIDENCE GRADING", first character = grade (1-5).
2. Position parsing: protein HGVS from the variant column (inhA_p.<...>). Missense /
   single-residue entries: score = min heavy-atom distance of that residue to ZID (1ZID).
   In-frame deletions/delins spanning a..b: min over the span. Insertions/dups between a and b:
   min over the two flanking residues (mirrors the rpob addendum convention).
3. Failure rule: an entry is a G1 failure if its residue (or entire span) is unresolved in 1ZID.
4. Swap-class operational definitions (for locked rule R3): formal charges at pH 7 taken as
   R=+1, K=+1, H=0, D=-1, E=-1, all others 0; "net charge change" = mutant charge != WT charge.
   Residue volumes (A^3, Zamyatnin): A 88.6, R 173.4, N 114.1, D 111.1, C 108.5, Q 143.8,
   E 138.4, G 60.1, H 153.2, I 166.7, L 166.7, K 168.6, M 162.9, F 189.9, P 112.7, S 89.0,
   T 116.1, W 227.8, Y 193.6, V 140.0; R3 triggers when RSA < 0.20 AND (charge change OR
   Pro gain/loss OR |delta volume| > 80 A^3). RSA from Shrake-Rupley on the protein-only 1ZID
   chain A (Biopython, default probe).
5. Effect types present in the WHO inhA subset (both drugs): missense_variant,
   synonymous_variant, upstream_gene_variant only. No stop_gained/frameshift/in-frame indel rows
   exist, so locked rule R0 has no instances; noted for completeness.
6. Structure-sequence fidelity (G2-fidelity result, observed while building the map): 1ZID
   chain A ATOM records span residues 2-269 and match P9WGR1 except position 2 (Ala in
   structure, Thr in reference) - an N-terminal cloning/engineering artifact consistent with the
   1ZID COMPND MUTATION flag. Position 2 is excluded from atlas claims. No other mismatches.
   WHO/H37Rv numbering == P9WGR1 numbering (offset 0), so no renumbering map is needed beyond
   this check.
7. G0 result observed while building the map (before G1): S94 3.48 A, Y158 3.80 A,
   K165 2.88 A to the ZID adduct - all < 6.0 A -> G0 PASS.
8. Resolution cross-check (supporting, not gated): per-residue ZID distances from 2IDZ (2.0 A
   WT InhA + ZID) correlate with 1ZID at r = 0.997, mean |delta| 0.40 A over 268 residues.
9. Promoter line (descriptive, not gated): upstream_gene_variant entries are parsed for c.-NN
   position; entries at -1..-50 relative to the inhA start codon are the locked promoter class
   (REGULATORY-R candidates). Their WHO grade distribution is reported; they never enter the G1
   structural denominator.
