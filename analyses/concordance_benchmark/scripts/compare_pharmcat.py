"""Compare afripharmagen star-allele calls against PharmCAT v3.4.0 output.

Parses PharmCAT JSON reports and afripharmagen benchmark JSONL, then produces
a per-allele carrier detection comparison table (Table 2 of the manuscript).

Usage:
    python compare_pharmcat.py \
        --pharmcat-dir /path/to/pharmcat_results/reports/ \
        --afripharmagen /path/to/afripharmagen_benchmark.jsonl \
        --samples /path/to/1000g_african_samples.tsv \
        --catalog /path/to/african_alleles.json \
        --output-dir /path/to/results/
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


def load_african_sample_ids(samples_path: Path) -> set[str]:
    """Load African sample IDs from the TSV."""
    ids = set()
    with open(samples_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ids.add(row["sample_id"])
    return ids


def load_catalog_alleles(catalog_path: Path) -> dict[str, list[str]]:
    """Load alleles to track from the catalog JSON."""
    with open(catalog_path) as f:
        data = json.load(f)
    alleles: dict[str, list[str]] = defaultdict(list)
    for entry in data["entries"]:
        alleles[entry["gene"]].append(entry["allele_name"])
    return dict(alleles)


def _init_stats(catalog_alleles: dict[str, list[str]]) -> dict[str, dict[str, dict]]:
    return {
        gene: {allele: {"carriers": 0, "no_call": 0, "total": 0} for allele in alleles}
        for gene, alleles in catalog_alleles.items()
    }


def _update_allele_stats(
    stats: dict[str, dict[str, dict]],
    gene: str,
    alleles: list[str],
    gene_data: dict,
) -> None:
    diplotypes = gene_data.get("sourceDiplotypes", [])
    if not diplotypes:
        for allele in alleles:
            stats[gene][allele]["no_call"] += 1
            stats[gene][allele]["total"] += 1
        return
    label = diplotypes[0].get("label", "Unknown/Unknown")
    for allele in alleles:
        stats[gene][allele]["total"] += 1
        if "Unknown" in label:
            stats[gene][allele]["no_call"] += 1
        elif allele in label:
            stats[gene][allele]["carriers"] += 1


def _process_report_file(
    report_file: Path,
    african_ids: set[str],
    catalog_alleles: dict[str, list[str]],
    stats: dict[str, dict[str, dict]],
) -> None:
    sample_id = report_file.name.split(".")[2]
    if sample_id not in african_ids:
        return
    with open(report_file) as f:
        data = json.load(f)
    for gene, alleles in catalog_alleles.items():
        gene_data = data.get("genes", {}).get(gene)
        if gene_data is not None:
            _update_allele_stats(stats, gene, alleles, gene_data)


def parse_reference_reports(  # noqa: SC200
    report_dir: Path, african_ids: set[str], catalog_alleles: dict[str, list[str]]
) -> dict[str, dict[str, dict]]:
    """Parse reference-caller JSON reports for African samples.

    Returns: {gene: {allele: {"carriers": int, "no_call": int, "total": int}}}
    """
    stats = _init_stats(catalog_alleles)
    for report_file in report_dir.glob("*.report.json"):
        _process_report_file(report_file, african_ids, catalog_alleles, stats)
    return stats


def parse_afripharmagen_results(
    jsonl_path: Path, catalog_alleles: dict[str, list[str]]
) -> dict[str, dict[str, int]]:
    """Parse afripharmagen JSONL for carrier counts.

    Returns: {gene: {allele: carrier_count}}
    """
    counts: dict[str, dict[str, int]] = {
        gene: dict.fromkeys(names, 0)
        for gene, names in catalog_alleles.items()
    }

    with open(jsonl_path) as f:
        for line in f:
            r = json.loads(line)
            if r.get("error"):
                continue
            gene = r["gene"]
            if gene not in catalog_alleles:
                continue
            diplotype = r.get("diplotype", "")
            for allele in catalog_alleles[gene]:
                if allele in diplotype:
                    counts[gene][allele] += 1

    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare PharmCAT vs afripharmagen")
    parser.add_argument("--pharmcat-dir", required=True, type=Path)
    parser.add_argument("--afripharmagen", required=True, type=Path)
    parser.add_argument(
        "--samples", type=Path,
        default=Path(__file__).resolve().parent.parent.parent.parent
        / "data" / "samples" / "1000g_african_samples.tsv",
    )
    parser.add_argument(
        "--catalog", type=Path,
        default=Path(__file__).resolve().parent.parent.parent.parent
        / "data" / "alleles" / "african_alleles.json",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()

    african_ids = load_african_sample_ids(args.samples)
    catalog_alleles = load_catalog_alleles(args.catalog)

    print(f"African samples: {len(african_ids)}")
    print(f"Catalog alleles: {sum(len(v) for v in catalog_alleles.values())} across {len(catalog_alleles)} genes")
    print()

    print("Parsing PharmCAT reports...")
    pc_stats = parse_reference_reports(args.pharmcat_dir, african_ids, catalog_alleles)

    print("Parsing afripharmagen results...")
    afri_counts = parse_afripharmagen_results(args.afripharmagen, catalog_alleles)

    # Output
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Print and write comparison table
    print()
    print("=" * 80)
    print("PHARMCAT v3.4.0 vs AFRIPHARMAGEN (n=661 African samples)")
    print("=" * 80)
    print()
    header = "Gene".ljust(10) + " " + "Allele".ljust(8) + " " + "afripharmagen".ljust(15) + " " + "PharmCAT".ljust(12) + " " + "PC_NoCalls".ljust(12) + " " + "Agreement"  # noqa: SC200
    print(header)
    print("-" * 75)

    report = {"n_samples": len(african_ids), "comparisons": []}

    for gene, alleles in catalog_alleles.items():
        for allele in alleles:
            afri = afri_counts[gene][allele]
            pc = pc_stats[gene][allele]["carriers"]
            no_call = pc_stats[gene][allele]["no_call"]

            if afri > 0 and pc > 0:
                agreement = f"{min(afri, pc) / max(afri, pc) * 100:.1f}%"
            elif afri == 0 and pc == 0:
                agreement = "N/A (none detected)"
            elif no_call == 661:
                agreement = "PharmCAT: NO CALL"
            else:
                agreement = f"0% ({pc} vs {afri})"

            print(f"{gene:<10} {allele:<8} {afri:<15} {pc:<12} {no_call:<12} {agreement}")

            report["comparisons"].append({
                "gene": gene,
                "allele": allele,
                "afripharmagen_carriers": afri,
                "pharmcat_carriers": pc,
                "pharmcat_no_calls": no_call,
                "agreement": agreement,
            })

    # Write JSON report
    json_path = args.output_dir / "pharmcat_comparison.json"
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)

    # Write TSV
    tsv_path = args.output_dir / "pharmcat_comparison.tsv"
    with open(tsv_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", "Allele", "afripharmagen_carriers", "PharmCAT_carriers", "PharmCAT_no_calls", "Agreement"])
        for c in report["comparisons"]:
            w.writerow([c["gene"], c["allele"], c["afripharmagen_carriers"],
                       c["pharmcat_carriers"], c["pharmcat_no_calls"], c["agreement"]])

    print()
    print(f"JSON: {json_path}")
    print(f"TSV:  {tsv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
