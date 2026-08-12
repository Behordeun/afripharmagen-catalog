"""Per-population allele frequency estimation from benchmark diplotype calls.

Computes carrier frequencies and allele frequencies for every allele in the
catalog across each 1000 Genomes African population. Identifies cases where:
  - A clinically actionable allele crosses the 1% clinical relevance threshold
    in some populations but not others
  - Aggregate "AFR" frequency masks population-specific gradients

Assumes diploid genomes. Allele frequency = (2*homozygotes + heterozygotes) / (2*N).
Carrier frequency = (homozygotes + heterozygotes) / N.

Outputs:
  - TSV: per-population frequency table (manuscript-ready)
  - JSON: structured report for platform consumption
  - Summary: gradient analysis (max-min spread, clinical threshold crossings)

Usage:
    python estimate_frequencies.py \
        --benchmark /path/to/full_benchmark.jsonl \
        --catalog /path/to/african_alleles.json \
        --samples /path/to/1000g_african_samples.tsv \
        --output-dir /path/to/results/
"""

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path


def wilson_ci(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if total == 0:
        return (0.0, 0.0)
    p_hat = successes / total
    denom = 1 + z**2 / total
    center = (p_hat + z**2 / (2 * total)) / denom
    spread = (
        z * math.sqrt((p_hat * (1 - p_hat) + z**2 / (4 * total)) / total) / denom
    )
    return (max(0.0, center - spread), min(1.0, center + spread))


def load_catalog_alleles(catalog_path: Path) -> dict[str, list[str]]:
    """Load tracked alleles from catalog. Returns gene -> [allele_names]."""
    with open(catalog_path, encoding="utf-8") as f:
        data = json.load(f)
    alleles: dict[str, list[str]] = defaultdict(list)
    for entry in data["entries"]:
        alleles[entry["gene"]].append(entry["allele_name"])
    return dict(alleles)


def load_sample_populations(samples_path: Path) -> dict[str, str]:
    """Load sample -> population mapping from TSV."""
    mapping = {}
    with open(samples_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["sample_id"]] = row["population"]
    return mapping


def load_benchmark(path: Path) -> list[dict]:
    """Load benchmark JSONL results."""
    results = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                results.append(json.loads(line))
    return results


def extract_alleles_from_diplotype(diplotype: str) -> list[str]:
    """Parse 'CYP2D6 *17/*2' into ['*17', '*2'] (preserving duplicates for homozygosity)."""
    if not diplotype:
        return []
    parts = diplotype.strip().split()
    allele_part = parts[1] if len(parts) == 2 else diplotype
    return [a.strip() for a in allele_part.split("/") if a.strip()]


@dataclass
class PopulationFrequency:
    """Frequency data for one allele in one population."""

    population: str
    n_samples: int = 0
    n_chromosomes: int = 0
    allele_count: int = 0
    carrier_count: int = 0
    homozygote_count: int = 0
    allele_frequency: float = 0.0
    carrier_frequency: float = 0.0
    ci_lower: float = 0.0
    ci_upper: float = 0.0


@dataclass
class AlleleFrequencyResult:
    """Frequency data for one allele across all populations."""

    gene: str
    allele: str
    total_carriers: int = 0
    total_samples: int = 0
    aggregate_frequency: float = 0.0
    max_frequency: float = 0.0
    min_frequency: float = 0.0
    max_population: str = ""
    min_population: str = ""
    gradient: float = 0.0
    crosses_1pct_threshold: bool = False
    populations: list = field(default_factory=list)


def _init_counts(
    catalog_alleles: dict[str, list[str]],
) -> dict[tuple[str, str], dict[str, dict[str, int]]]:
    counts: dict[tuple[str, str], dict[str, dict[str, int]]] = {}
    for gene, alleles in catalog_alleles.items():
        for allele in alleles:
            counts[(gene, allele)] = defaultdict(
                lambda: {"allele_count": 0, "carrier_count": 0, "hom": 0}
            )
    return counts


def _process_freq_result(
    result: dict,
    catalog_alleles: dict[str, list[str]],
    sample_populations: dict[str, str],
    counts: dict[tuple[str, str], dict[str, dict[str, int]]],
    seen: dict[tuple[str, str], set[str]],
) -> None:
    if result.get("error"):
        return
    gene = result["gene"].upper()
    if gene not in catalog_alleles:
        return
    sample_id = result["sample_id"]
    population = result.get("population") or sample_populations.get(sample_id, "UNK")
    if sample_id in seen[(gene, population)]:
        return
    seen[(gene, population)].add(sample_id)
    diplotype_alleles = extract_alleles_from_diplotype(result.get("diplotype", ""))
    for allele in catalog_alleles[gene]:
        n = diplotype_alleles.count(allele)
        if n > 0:
            pop_counts = counts[(gene, allele)][population]
            pop_counts["allele_count"] += n
            pop_counts["carrier_count"] += 1
            if n == 2:
                pop_counts["hom"] += 1


def _build_pop_frequency(
    pop: str,
    n_samples: int,
    pop_data: dict[str, int],
) -> PopulationFrequency:
    n_chromosomes = 2 * n_samples
    af = pop_data["allele_count"] / n_chromosomes if n_chromosomes > 0 else 0
    cf = pop_data["carrier_count"] / n_samples if n_samples > 0 else 0
    ci_low, ci_high = wilson_ci(pop_data["allele_count"], n_chromosomes)
    return PopulationFrequency(
        population=pop,
        n_samples=n_samples,
        n_chromosomes=n_chromosomes,
        allele_count=pop_data["allele_count"],
        carrier_count=pop_data["carrier_count"],
        homozygote_count=pop_data["hom"],
        allele_frequency=af,
        carrier_frequency=cf,
        ci_lower=ci_low,
        ci_upper=ci_high,
    )


def _build_allele_result(
    gene: str,
    allele: str,
    pop_sizes: dict[str, int],
    counts: dict[tuple[str, str], dict[str, dict[str, int]]],
) -> AlleleFrequencyResult:
    ar = AlleleFrequencyResult(gene=gene, allele=allele)
    empty: dict[str, int] = {"allele_count": 0, "carrier_count": 0, "hom": 0}
    pop_freqs = [
        _build_pop_frequency(pop, pop_sizes[pop], counts[(gene, allele)].get(pop, empty))
        for pop in sorted(pop_sizes)
    ]
    ar.total_carriers = sum(pf.carrier_count for pf in pop_freqs)
    ar.total_samples = sum(pop_sizes.values())
    ar.populations = pop_freqs
    total_allele_count = sum(pf.allele_count for pf in pop_freqs)
    total_chromosomes = sum(pf.n_chromosomes for pf in pop_freqs)
    ar.aggregate_frequency = total_allele_count / total_chromosomes if total_chromosomes > 0 else 0
    if pop_freqs:
        max_freq, max_pop = max((pf.allele_frequency, pf.population) for pf in pop_freqs)
        min_freq, min_pop = min((pf.allele_frequency, pf.population) for pf in pop_freqs)
        ar.max_frequency, ar.max_population = max_freq, max_pop
        ar.min_frequency, ar.min_population = min_freq, min_pop
        ar.gradient = max_freq - min_freq
        ar.crosses_1pct_threshold = (
            any(pf.allele_frequency >= 0.01 for pf in pop_freqs)
            and any(pf.allele_frequency < 0.01 for pf in pop_freqs)
        )
    return ar


def compute_frequencies(
    benchmark: list[dict],
    catalog_alleles: dict[str, list[str]],
    sample_populations: dict[str, str],
) -> list[AlleleFrequencyResult]:
    """Compute per-population allele frequencies from benchmark diplotype calls."""
    pop_sizes: dict[str, int] = defaultdict(int)
    for pop in sample_populations.values():
        pop_sizes[pop] += 1

    counts = _init_counts(catalog_alleles)
    seen: dict[tuple[str, str], set[str]] = defaultdict(set)

    for result in benchmark:
        _process_freq_result(result, catalog_alleles, sample_populations, counts, seen)

    return [
        _build_allele_result(gene, allele, pop_sizes, counts)
        for gene, alleles in catalog_alleles.items()
        for allele in alleles
    ]


def write_tsv(results: list[AlleleFrequencyResult], output_path: Path) -> None:
    """Write per-population frequency table."""
    # Collect all populations
    all_pops = sorted({pf.population for ar in results for pf in ar.populations})

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")

        # Header
        header = ["Gene", "Allele", "Aggregate_AF"]
        for pop in all_pops:
            header.extend([f"{pop}_AF", f"{pop}_CI_Low", f"{pop}_CI_High", f"{pop}_N"])
        header.extend(["Max_Pop", "Min_Pop", "Gradient", "Crosses_1pct"])
        w.writerow(header)

        for ar in results:
            row = [ar.gene, ar.allele, f"{ar.aggregate_frequency:.4f}"]
            pop_lookup = {pf.population: pf for pf in ar.populations}
            for pop in all_pops:
                pf = pop_lookup.get(pop)
                if pf:
                    row.extend([
                        f"{pf.allele_frequency:.4f}",
                        f"{pf.ci_lower:.4f}",
                        f"{pf.ci_upper:.4f}",
                        str(pf.n_samples),
                    ])
                else:
                    row.extend(["0.0000", "0.0000", "0.0000", "0"])
            row.extend([
                ar.max_population, ar.min_population,
                f"{ar.gradient:.4f}", str(ar.crosses_1pct_threshold),
            ])
            w.writerow(row)


def write_json(results: list[AlleleFrequencyResult], output_path: Path) -> None:
    """Write structured JSON report."""
    report = {
        "catalog_version": "0.1.0",
        "total_populations": len({pf.population for ar in results for pf in ar.populations}),
        "total_alleles": len(results),
        "alleles": [],
    }
    for ar in results:
        entry = {
            "gene": ar.gene,
            "allele": ar.allele,
            "aggregate_frequency": ar.aggregate_frequency,
            "total_carriers": ar.total_carriers,
            "total_samples": ar.total_samples,
            "gradient": ar.gradient,
            "max_frequency": ar.max_frequency,
            "max_population": ar.max_population,
            "min_frequency": ar.min_frequency,
            "min_population": ar.min_population,
            "crosses_1pct_threshold": ar.crosses_1pct_threshold,
            "populations": [asdict(pf) for pf in ar.populations],
        }
        report["alleles"].append(entry)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


def print_summary(results: list[AlleleFrequencyResult]) -> None:
    """Print gradient analysis summary."""
    print(f"\n{'='*80}")
    print("POPULATION FREQUENCY ANALYSIS")
    print(f"{'='*80}\n")

    print(f"{'Gene':<10} {'Allele':<8} {'Aggregate':<11} {'Max (Pop)':<16} "
          f"{'Min (Pop)':<16} {'Gradient':<10} {'1% Cross'}")
    print("-" * 80)

    for ar in results:
        max_str = f"{ar.max_frequency:.3f} ({ar.max_population})"
        min_str = f"{ar.min_frequency:.3f} ({ar.min_population})"
        cross_str = "YES" if ar.crosses_1pct_threshold else "no"
        print(f"{ar.gene:<10} {ar.allele:<8} {ar.aggregate_frequency:<11.4f} "
              f"{max_str:<16} {min_str:<16} {ar.gradient:<10.4f} {cross_str}")

    # Highlight gradients > 10 percentage points
    large_gradients = [ar for ar in results if ar.gradient >= 0.10]
    if large_gradients:
        print(f"\n{'─'*60}")
        print("CLINICALLY SIGNIFICANT GRADIENTS (>10 percentage points):")
        print(f"{'─'*60}")
        for ar in large_gradients:
            print(f"  {ar.gene} {ar.allele}: {ar.min_frequency:.1%} ({ar.min_population}) "
                  f"to {ar.max_frequency:.1%} ({ar.max_population}) "
                  f"[spread: {ar.gradient:.1%}]")

    # Threshold crossings
    crossings = [ar for ar in results if ar.crosses_1pct_threshold]
    if crossings:
        print(f"\n{'─'*60}")
        print("ALLELES CROSSING 1% CLINICAL THRESHOLD BETWEEN POPULATIONS:")
        print(f"{'─'*60}")
        for ar in crossings:
            above = [pf.population for pf in ar.populations if pf.allele_frequency >= 0.01]
            below = [pf.population for pf in ar.populations if pf.allele_frequency < 0.01]
            print(f"  {ar.gene} {ar.allele}: above 1% in {', '.join(above)}; "
                  f"below 1% in {', '.join(below)}")

    print()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Per-population allele frequency estimation"
    )
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument(
        "--catalog", type=Path,
        default=Path(__file__).resolve().parent.parent.parent.parent
        / "data" / "alleles" / "african_alleles.json",
    )
    parser.add_argument(
        "--samples", type=Path,
        default=Path(__file__).resolve().parent.parent.parent.parent
        / "data" / "samples" / "1000g_african_samples.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()

    for path, name in [(args.benchmark, "benchmark"), (args.catalog, "catalog"),
                       (args.samples, "samples")]:
        if not path.exists():
            print(f"ERROR: {name} file not found: {path}", file=sys.stderr)
            return 1

    catalog_alleles = load_catalog_alleles(args.catalog)
    sample_populations = load_sample_populations(args.samples)
    benchmark = load_benchmark(args.benchmark)

    if not benchmark:
        print("ERROR: No results in benchmark file", file=sys.stderr)
        return 1

    results = compute_frequencies(benchmark, catalog_alleles, sample_populations)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.output_dir / "population_frequencies.tsv"
    json_path = args.output_dir / "population_frequencies.json"

    write_tsv(results, tsv_path)
    write_json(results, json_path)
    print_summary(results)

    print(f"TSV:  {tsv_path}")
    print(f"JSON: {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
