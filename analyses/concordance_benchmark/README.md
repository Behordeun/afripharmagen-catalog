# Concordance Benchmark

Evaluation of star-allele calling concordance between the afripharmagen pipeline and PharmCAT v2.13.0 across 661 African genomes from the 1000 Genomes Project.

## Method

1. Ground truth diplotypes assigned via direct variant inspection at catalog-defined positions using phased haplotype data from 1000 Genomes Project statistical phasing.
2. afripharmagen StarCaller run on each sample's VCF for all catalog genes.
3. Per-allele concordance computed as fraction of carriers correctly identified.

## Populations

| Code | Description | N |
|------|-------------|---|
| YRI | Yoruba in Ibadan, Nigeria | 108 |
| LWK | Luhya in Webuye, Kenya | 99 |
| GWD | Gambian in Western Divisions | 113 |
| MSL | Mende in Sierra Leone | 85 |
| ESN | Esan in Nigeria | 99 |
| ACB | African Caribbeans in Barbados | 96 |
| ASW | Americans of African Ancestry in SW USA | 61 |

## Results

Results will be deposited in `results/` upon benchmark completion.

## Scripts

- `scripts/run_benchmark.py` - Full benchmark execution
- `scripts/compute_concordance.py` - Per-allele concordance statistics
