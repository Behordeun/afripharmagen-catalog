# Allele Catalog Analysis

Construction, validation, and documentation of the African pharmacogenomic allele catalog.

## Catalog Schema

Each entry in `data/alleles/african_alleles.json` contains:

| Field | Type | Description |
|-------|------|-------------|
| gene | string | Pharmacogene name (e.g., CYP2D6) |
| allele_name | string | Star allele designation (e.g., *17) |
| defining_variants | array | rsID, GRCh38 position, ref, alt for each defining variant |
| function | string | Functional classification (decreased_function, no_function, normal_function) |
| activity_score | float | CPIC activity score (0.0, 0.5, or 1.0) |
| evidence_level | string | L1 (in vitro confirmed) or L2 (computationally predicted) |
| populations | array | 1000 Genomes population codes where allele is reported |
| frequency_range | string | Published frequency range across African populations |
| source_pmids | array | PubMed identifiers for source literature |
| in_pharmvar | boolean | Whether allele is cataloged in PharmVar |
| notes | string | Clinical context and curation notes |

## Current Catalog (v0.1.0)

8 alleles across 5 genes: CYP2D6*17, *29, *45; CYP2B6*6, *18; CYP2C19*9; CYP3A5*1; NAT2*14.

## Scripts

- `scripts/validate_catalog.py` - Schema validation and cross-reference checks
