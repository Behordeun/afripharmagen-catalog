# Concordance Benchmark

Comparison of afripharmagen star-allele calls against PharmCAT v3.4.0 on 661 African samples from the 1000 Genomes Project.

## Method

1. PharmCAT v3.4.0 run on all 3,202 1000 Genomes samples using the official VCF preprocessor and default configuration.
2. afripharmagen StarCaller run on the same VCFs for all 661 African samples across 10 genes.
3. Per-allele carrier detection compared between tools for all 9 catalog alleles.

## Key Finding

PharmCAT v3.4.0 cannot produce any CYP2D6 diplotype call from standard variant-only WGS VCFs (returns "Unknown/Unknown" for all 661 African samples). This is a documented limitation of PharmCAT's architecture: it requires explicit genotype calls at monomorphic reference positions that standard variant-calling pipelines do not emit. No prior study had quantified this limitation in terms of African-specific allele carrier counts and clinical drug-prescribing consequences.

afripharmagen's reduced-position matching strategy identifies 243 CYP2D6*17 carriers and 134 CYP2D6*29 carriers from the same input by requiring only the 1-2 core defining variants per allele (which are polymorphic and therefore present in the VCF).

For other genes (CYP2B6, CYP2C9, CYP2C19, NAT2), PharmCAT and afripharmagen show 95-100% concordance. For CYP3A5, PharmCAT outperforms afripharmagen by correctly identifying 65 CYP3A5*7 non-expressers that the current catalog misses.

## Scripts

- `scripts/compare_pharmcat.py` - Parses PharmCAT JSON reports and afripharmagen JSONL, produces Table 2

## Results

- `results/pharmcat_comparison.json` - Machine-readable comparison
- `results/pharmcat_comparison.tsv` - TSV for manuscript table

## Reproducing

```bash
python scripts/compare_pharmcat.py \
    --pharmcat-dir /path/to/pharmcat_results/reports/ \
    --afripharmagen ../../results/afripharmagen_benchmark.jsonl \
    --output-dir results/
```

## PharmCAT Version

- PharmCAT v3.4.0 (downloaded from https://github.com/PharmGKB/PharmCAT/releases/tag/v3.4.0)
- Run date: August 2026
