# afripharmagen-catalog

Curated pharmacogenomic data and reproducible analyses for sub-Saharan African populations.

## Contents

### data/

Core datasets shared across analyses:

- `alleles/african_alleles.json` - Machine-readable catalog of clinically actionable pharmacogenomic alleles with African-specific frequency ranges, functional annotations, and evidence levels.
- `samples/1000g_african_samples.tsv` - Population metadata for 661 samples from seven 1000 Genomes Project African populations (YRI, LWK, GWD, MSL, ESN, ACB, ASW).

### analyses/

Reproducible analysis modules. Each directory contains scripts and results for a specific line of investigation:

- `allele_catalog/` - Catalog construction, validation, and schema documentation
- `concordance_benchmark/` - Star-allele calling concordance evaluation against PharmCAT
- `population_frequencies/` - Per-population allele frequency estimation and gradients

## Data Sources

- Allele definitions derived from PharmVar, CPIC, and published literature (PMIDs documented per entry)
- Whole-genome sequence data from the 1000 Genomes Project high-coverage dataset (Byrska-Bishop et al., Cell, 2022)
- VCFs aligned to GRCh38, accessed via the International Genome Sample Resource (https://www.internationalgenome.org/)

## Citation

If you use this data or these analyses, please cite:

> Sulaiman MA, Oyeyemi BF. A Curated Pharmacogenomic Allele Catalog for Sub-Saharan African Populations. [Manuscript in preparation]. 2026.

## License

Data and analysis scripts are released under the MIT License. See `LICENSE` for details.
