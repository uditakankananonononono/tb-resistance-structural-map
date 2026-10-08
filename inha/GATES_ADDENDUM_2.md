# Gates addendum 2 - promoter-class mapping correction (descriptive line only; G1 untouched)
Written 2026-10-08 ~19:10 IST, after the grade-1/2 variant list was printed by the first
02_validate run (disclosed: this correction IS informed by observed outcome positions; it touches
only the descriptive promoter line, never the gated structural denominator or thresholds).

Addendum 1 item 9 operationally mis-mapped the locked promoter class. GATES_LOCKED.md defines
the class as "variants upstream of the fabG1 start codon (the fabG1-inhA operon promoter),
position class -1..-50 relative to the fabG1 start" - fabG1-relative. Addendum 1 mistakenly
implemented it as -1..-50 relative to the INHA start. In WHO catalogue c. numbering (relative to
the inhA start codon), the fabG1 start maps as follows, derived from the WHO master file's own
genomic position column: inhA c.1 is at genomic 1674202 (plus strand); the grade-1 entry
inhA_c.-777C>T (genomic 1673425) is the canonical fabG1 c.-15C>T operon-promoter mutation, so
fabG1 c.1 == inhA c.-762 and the locked -1..-50 fabG1-relative class == inhA-relative
c.-762..c.-811 (genomic 1673341..1673390 + offsets). The grade-1/2 promoter entries
(c.-770T>A/C/G, c.-777C>T, c.-778A>G, c.-779G>T) all fall inside this corrected class; the
addendum-1 class (-1..-50 inhA-relative) contains none of them and stays in the record as
written. The corrected descriptive line reports the fabG1-relative class. inhA_c.-154G>A (grade
1 INH / grade 2 ETH) is NOT in the operon-promoter class; it lies inside the fabG1 coding region
by this mapping and is reported separately as a descriptive out-of-class observation.
G1 (coding, structural) is unchanged: 1 unique grade-1/2 protein-altering variant (S94A).
