# Gates addendum 1 - operational definitions (written before G1 was computed)
Written 2026-09-23 ~13:45 IST. At this point only grade counts and variant names had been viewed; no distances had been joined to WHO entries.
Thresholds and the gate are unchanged. This only states how entries are scored:
1. Position parsing: HGVS protein positions in WHO numbering (1172-aa RpoB). Structure numbering was mapped by global alignment to UniProt P9WGY9 residues 7-1178 (100% identity over 1126 resolved residues; 5UHC auth numbering = WHO + 6).
2. Missense / single-residue entries: score = min heavy-atom distance of that residue to RFP.
3. In-frame deletions/delins spanning residues a..b: score = min over residues a..b. Insertions/dups between a and b: score = min over the two flanking residues (dup: the duplicated residue(s)).
4. An entry is a failure if its residue (or all residues of its span) is unresolved in 5UHC.
5. G0 result observed while building the map (before G1): S450 2.43 A, H445 3.42 A, D435 3.14 A -> G0 PASS.
