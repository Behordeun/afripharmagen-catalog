# Population Frequency Analysis

Per-population allele frequency estimation for catalog alleles across seven African populations.

## Method

- Single-variant alleles: frequency = proportion of chromosomes carrying the alternate allele.
- Multi-variant haplotype alleles (CYP2D6*17, *29): frequency estimated from statistically phased haplotypes.
- CYP3A5*1: estimated indirectly as 1 - frequency(CYP3A5*3, rs776746).
- Hardy-Weinberg equilibrium assumed for phenotype frequency derivation from allele frequencies.

## Scripts

- `scripts/estimate_frequencies.py` - Per-population allele frequency computation
- `scripts/plot_gradients.py` - Visualization of frequency gradients across populations
