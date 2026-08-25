# afripharmagen-catalog

Curated pharmacogenomic data and reproducible analyses for sub-Saharan African populations.

## Contents

### data/

Core datasets:

- `alleles/african_alleles.json` - Machine-readable catalog of 10 clinically actionable pharmacogenomic alleles across 6 genes, with African-specific frequency ranges, functional annotations, and evidence levels.
- `samples/1000g_african_samples.tsv` - Population metadata for 661 samples from seven 1000 Genomes Project African populations (YRI, LWK, GWD, MSL, ESN, ACB, ASW).

### results/

Benchmark outputs:

- `afripharmagen_benchmark.jsonl` - Star-allele calls from afripharmagen for all 661 samples x 10 genes (6,610 results).

### analyses/

Reproducible analysis modules:

- `concordance_benchmark/` - Comparison of afripharmagen vs. PharmCAT v3.4.0 carrier detection
- `population_frequencies/` - Per-population allele frequency estimation, gradient analysis, and Wilson score confidence intervals
- `allele_catalog/` - Catalog schema documentation

## Key Findings

PharmCAT v3.4.0 cannot call CYP2D6 from standard variant-only VCF output (0/661 African samples callable), a known VCF format limitation documented by the PharmCAT team. afripharmagen's reduced-position strategy identifies 243 CYP2D6\*17 and 134 CYP2D6\*29 carriers from the same input without reprocessing. For other genes, both tools agree (95-100%).

Three alleles show inter-population frequency gradients exceeding 10 percentage points across the seven African populations: CYP3A5\*1 (19.3pp), CYP2D6\*17 (13.5pp), and CYP2D6\*29 (13.1pp). Per-population frequency data with 95% Wilson confidence intervals are available in `analyses/population_frequencies/results/`.

## Data Sources

- Allele definitions derived from PharmVar v6.0.4 and published literature
- Whole-genome sequence data from the 1000 Genomes Project high-coverage dataset (Byrska-Bishop et al., Cell, 2022)
- PharmCAT v3.4.0 (https://github.com/PharmGKB/PharmCAT/releases/tag/v3.4.0)

## Quickstart

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/Behordeun/afripharmagen-catalog.git
cd afripharmagen-catalog

# Create environment and install dependencies
uv sync

# Run the full analysis pipeline
uv run python run_analyses.py
```

That's it. `uv sync` reads `pyproject.toml` and installs `numpy`, `scipy`, and `matplotlib` into a local `.venv`.

## Running All Analyses

The top-level pipeline script orchestrates the full workflow:

```bash
# Run everything except the PharmCAT concordance (requires external reports)
uv run python run_analyses.py

# Include concordance benchmark if you have PharmCAT reports
uv run python run_analyses.py --pharmcat-dir /path/to/pharmcat_results/reports/

# Skip figure generation (headless CI environments)
uv run python run_analyses.py --skip-plots
```

The pipeline runs in order: catalog validation, frequency estimation, HWE testing, gradient plots, and (optionally) PharmCAT concordance.

## Dependencies

All dependencies are declared in `pyproject.toml` and pinned to compatible ranges:

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | >=1.24, <3.0 | Array operations, frequency calculations |
| scipy | >=1.11, <2.0 | Chi-squared HWE tests |
| matplotlib | >=3.8, <4.0 | Frequency gradient figures |

No additional packages are needed to reproduce all analyses.

## Citation

See [`CITATION.cff`](CITATION.cff) for citation information.

## License

MIT License. See `LICENSE`.
