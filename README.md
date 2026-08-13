# afripharmagen-catalog

Curated pharmacogenomic data and reproducible analyses for sub-Saharan African populations.

## Contents

### data/

Core datasets:

- `alleles/african_alleles.json` - Machine-readable catalog of 9 clinically actionable pharmacogenomic alleles across 6 genes, with African-specific frequency ranges, functional annotations, and evidence levels.
- `samples/1000g_african_samples.tsv` - Population metadata for 661 samples from seven 1000 Genomes Project African populations (YRI, LWK, GWD, MSL, ESN, ACB, ASW).

### results/

Benchmark outputs:

- `afripharmagen_benchmark.jsonl` - Star-allele calls from afripharmagen for all 661 samples x 10 genes (6,610 results).

### analyses/

Reproducible analysis modules:

- `concordance_benchmark/` - Comparison of afripharmagen vs. PharmCAT v3.4.0 carrier detection
- `population_frequencies/` - Per-population allele frequency estimation and gradient analysis
- `allele_catalog/` - Catalog schema documentation

## Key Finding

PharmCAT v3.4.0 cannot call CYP2D6 from standard variant-only VCF output (0/661 African samples callable), a known VCF format limitation documented by the PharmCAT team. afripharmagen's reduced-position strategy identifies 243 CYP2D6\*17 and 134 CYP2D6\*29 carriers from the same input without reprocessing. For other genes, both tools agree (95-100%). For CYP3A5, PharmCAT correctly handles \*7 non-expressers that the current afripharmagen catalog misses.

## Data Sources

- Allele definitions derived from PharmVar v6.0.4 and published literature
- Whole-genome sequence data from the 1000 Genomes Project high-coverage dataset (Byrska-Bishop et al., Cell, 2022)
- PharmCAT v3.4.0 (https://github.com/PharmGKB/PharmCAT/releases/tag/v3.4.0)

## Citation

See [`CITATION.cff`](CITATION.cff) for citation information.

## License

MIT License. See `LICENSE`.
