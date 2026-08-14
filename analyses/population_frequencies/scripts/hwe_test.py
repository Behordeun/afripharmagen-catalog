"""
Hardy-Weinberg Equilibrium test for per-population allele frequencies.

Tests each allele-population combination using a chi-squared goodness-of-fit
test (1 df) comparing observed genotype counts against HWE expectations.

Input: afripharmagen benchmark JSONL (diplotype calls)
Output: Summary of HWE test results with Bonferroni correction for 63 tests.

Usage:
    python analyses/population_frequencies/scripts/hwe_test.py
"""

import json
import re
from collections import defaultdict
from pathlib import Path

from scipy import stats

BENCHMARK_JSONL = Path(__file__).resolve().parent.parent.parent.parent / "results" / "afripharmagen_benchmark.jsonl"

TARGET_ALLELES = [
    ("CYP3A5", "*1"),
    ("CYP3A5", "*7"),
    ("CYP2D6", "*17"),
    ("CYP2D6", "*29"),
    ("CYP2B6", "*6"),
    ("CYP2B6", "*18"),
    ("CYP2C9", "*8"),
    ("CYP2C19", "*9"),
    ("NAT2", "*14"),
]

POPULATIONS = ["YRI", "ESN", "GWD", "MSL", "LWK", "ACB", "ASW"]


def count_genotypes(benchmark_path: Path) -> dict:
    """Count homozygous, heterozygous, and reference genotypes per allele per population."""
    genotype_counts = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: {"hom": 0, "het": 0, "ref": 0}))
    )

    with open(benchmark_path) as f:
        for line in f:
            r = json.loads(line)
            gene = r["gene"]
            pop = r["population"]
            diplotype = r["diplotype"]

            match = re.match(r"(\w+)\s+\*(.+)/\*(.+)", diplotype)
            if not match:
                continue
            a1, a2 = match.group(2), match.group(3)

            for tg, ta in TARGET_ALLELES:
                if gene == tg:
                    target = ta.lstrip("*")
                    copies = (1 if a1 == target else 0) + (1 if a2 == target else 0)
                    if copies == 2:
                        genotype_counts[tg][ta][pop]["hom"] += 1
                    elif copies == 1:
                        genotype_counts[tg][ta][pop]["het"] += 1
                    else:
                        genotype_counts[tg][ta][pop]["ref"] += 1

    return genotype_counts


def run_hwe_tests(genotype_counts: dict) -> list[dict]:
    """Run chi-squared HWE test for each allele-population pair."""
    results = []

    for gene, allele in TARGET_ALLELES:
        for pop in POPULATIONS:
            counts = genotype_counts[gene][allele][pop]
            n = counts["hom"] + counts["het"] + counts["ref"]
            if n == 0:
                continue

            p = (2 * counts["hom"] + counts["het"]) / (2 * n)
            q = 1 - p

            exp_hom = p**2 * n
            exp_het = 2 * p * q * n
            exp_ref = q**2 * n

            # Skip test if expected counts are too low for chi-squared
            if exp_hom < 0.5 and exp_het < 0.5:
                chi2 = 0.0
                pval = 1.0
            else:
                obs = [counts["hom"], counts["het"], counts["ref"]]
                exp = [exp_hom, exp_het, exp_ref]
                chi2 = sum((o - e) ** 2 / e for o, e in zip(obs, exp) if e > 0)
                pval = 1 - stats.chi2.cdf(chi2, df=1)

            results.append({
                "gene": gene,
                "allele": allele,
                "population": pop,
                "n": n,
                "obs_hom": counts["hom"],
                "obs_het": counts["het"],
                "obs_ref": counts["ref"],
                "allele_freq": p,
                "exp_hom": exp_hom,
                "exp_het": exp_het,
                "exp_ref": exp_ref,
                "chi2": chi2,
                "pvalue": pval,
            })

    return results


def main():
    if not BENCHMARK_JSONL.exists():
        print(f"ERROR: Benchmark file not found: {BENCHMARK_JSONL}")
        return 1

    print(f"Input: {BENCHMARK_JSONL}")
    print()

    genotype_counts = count_genotypes(BENCHMARK_JSONL)
    results = run_hwe_tests(genotype_counts)

    total_tests = len(results)
    bonferroni_threshold = 0.05 / total_tests
    pvalues = [r["pvalue"] for r in results]

    nominal = [r for r in results if r["pvalue"] < 0.05]
    significant = [r for r in results if r["pvalue"] < bonferroni_threshold]

    print(f"Total allele-population tests: {total_tests}")
    print(f"Bonferroni-corrected threshold: p < {bonferroni_threshold:.6f}")
    print(f"Significant after Bonferroni: {len(significant)}")
    print(f"Nominal p < 0.05 (uncorrected): {len(nominal)}")
    print(f"Minimum p-value: {min(pvalues):.6f}")
    print()

    if nominal:
        print("Tests with nominal p < 0.05:")
        for r in sorted(nominal, key=lambda x: x["pvalue"]):
            bf_flag = " [would not survive Bonferroni]"
            if r["pvalue"] < bonferroni_threshold:
                bf_flag = " *** SIGNIFICANT after Bonferroni ***"
            print(
                f"  {r['gene']} {r['allele']} {r['population']}: "
                f"n={r['n']}, genotypes=[{r['obs_hom']}hom, {r['obs_het']}het, {r['obs_ref']}ref], "
                f"p={r['allele_freq']:.3f}, chi2={r['chi2']:.3f}, "
                f"pval={r['pvalue']:.6f}{bf_flag}"
            )
    else:
        print("No tests reached nominal p < 0.05.")

    print()
    print("CONCLUSION: ", end="")
    if len(significant) == 0:
        print(
            f"No significant HWE deviations after Bonferroni correction "
            f"({len(nominal)}/{total_tests} tests reached nominal p<0.05; "
            f"none survived correction at p<{bonferroni_threshold:.4f})."
        )
    else:
        print(
            f"{len(significant)} significant HWE deviations detected after Bonferroni correction."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
