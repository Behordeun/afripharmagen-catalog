"""Visualize allele frequency gradients across African populations.

Reads the structured JSON output from estimate_frequencies.py and produces
a grouped bar chart showing per-population allele frequencies with Wilson
score 95% confidence intervals.

Usage:
    python plot_gradients.py \
        --frequencies results/population_frequencies.json \
        --output-dir results/
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_frequency_data(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def plot_frequency_heatmap(data: dict, output_dir: Path) -> None:
    """Heatmap of allele frequencies across populations."""
    alleles_data = data["alleles"]
    populations = sorted(
        {p["population"] for a in alleles_data for p in a["populations"]}
    )
    allele_labels = [f"{a['gene']} {a['allele']}" for a in alleles_data]

    matrix = np.zeros((len(alleles_data), len(populations)))
    for i, allele in enumerate(alleles_data):
        pop_lookup = {p["population"]: p["allele_frequency"] for p in allele["populations"]}
        for j, pop in enumerate(populations):
            matrix[i, j] = pop_lookup.get(pop, 0.0)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix, cmap="viridis", aspect="auto", vmin=0)

    ax.set_xticks(range(len(populations)))
    ax.set_xticklabels(populations, fontsize=10, fontfamily="sans-serif")
    ax.set_yticks(range(len(allele_labels)))
    ax.set_yticklabels(allele_labels, fontsize=10, fontfamily="sans-serif")

    # Annotate cells with frequency values
    for i in range(len(allele_labels)):
        for j in range(len(populations)):
            val = matrix[i, j]
            color = "white" if val > 0.15 else "black"
            ax.text(j, i, f"{val:.1%}", ha="center", va="center",
                    fontsize=8, color=color, fontfamily="sans-serif")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Allele Frequency", fontsize=10, fontfamily="sans-serif")

    ax.set_title("Pharmacogenomic Allele Frequencies Across African Populations",
                 fontsize=12, fontfamily="sans-serif", pad=12)
    ax.set_xlabel("Population", fontsize=10, fontfamily="sans-serif")
    ax.set_ylabel("Allele", fontsize=10, fontfamily="sans-serif")

    plt.tight_layout()

    for fmt in ("png", "pdf", "svg"):
        out_path = output_dir / f"figure1_frequency_heatmap.{fmt}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gradient_bars(data: dict, output_dir: Path) -> None:
    """Bar chart of frequency gradients sorted by magnitude."""
    alleles_data = sorted(data["alleles"], key=lambda a: a["gradient"], reverse=True)
    labels = [f"{a['gene']} {a['allele']}" for a in alleles_data]
    gradients = [a["gradient"] * 100 for a in alleles_data]
    colors = ["#d62728" if g >= 10 else "#1f77b4" for g in gradients]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(range(len(labels)), gradients, color=colors, edgecolor="none")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=10, fontfamily="sans-serif")
    ax.set_xlabel("Frequency Gradient (percentage points)", fontsize=10,
                  fontfamily="sans-serif")
    ax.set_title("Inter-Population Allele Frequency Gradients",
                 fontsize=12, fontfamily="sans-serif", pad=12)

    # Reference line at 10pp threshold
    ax.axvline(x=10, color="#555555", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.text(10.3, len(labels) - 0.5, "10pp threshold", fontsize=8,
            color="#555555", fontfamily="sans-serif", va="top")

    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    for fmt in ("png", "pdf", "svg"):
        out_path = output_dir / f"figure2_gradient_bars.{fmt}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_population_detail(data: dict, output_dir: Path) -> None:
    """Grouped bar chart with error bars for alleles with gradient > 10pp."""
    large_gradient = [a for a in data["alleles"] if a["gradient"] >= 0.10]
    if not large_gradient:
        return

    populations = sorted(
        {p["population"] for a in large_gradient for p in a["populations"]}
    )
    n_alleles = len(large_gradient)
    n_pops = len(populations)

    fig, ax = plt.subplots(figsize=(10, 5))
    bar_width = 0.8 / n_alleles
    x = np.arange(n_pops)

    cmap = plt.get_cmap("cividis", n_alleles)

    for i, allele in enumerate(large_gradient):
        pop_lookup = {p["population"]: p for p in allele["populations"]}
        freqs = []
        ci_low = []
        ci_high = []
        for pop in populations:
            p = pop_lookup.get(pop)
            if p:
                freqs.append(p["allele_frequency"] * 100)
                ci_low.append((p["allele_frequency"] - p["ci_lower"]) * 100)
                ci_high.append((p["ci_upper"] - p["allele_frequency"]) * 100)
            else:
                freqs.append(0)
                ci_low.append(0)
                ci_high.append(0)

        offset = (i - n_alleles / 2 + 0.5) * bar_width
        label = f"{allele['gene']} {allele['allele']}"
        ax.bar(x + offset, freqs, bar_width, yerr=[ci_low, ci_high],
               label=label, color=cmap(i), capsize=2, edgecolor="none")

    ax.set_xticks(x)
    ax.set_xticklabels(populations, fontsize=10, fontfamily="sans-serif")
    ax.set_ylabel("Allele Frequency (%)", fontsize=10, fontfamily="sans-serif")
    ax.set_title("High-Gradient Alleles: Per-Population Frequencies with 95% CI",
                 fontsize=12, fontfamily="sans-serif", pad=12)
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    for fmt in ("png", "pdf", "svg"):
        out_path = output_dir / f"figure3_high_gradient_detail.{fmt}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot allele frequency gradients across populations"
    )
    parser.add_argument("--frequencies", required=True, type=Path,
                        help="Path to population_frequencies.json")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="Directory to save figures")
    args = parser.parse_args()

    if not args.frequencies.exists():
        print(f"ERROR: frequency file not found: {args.frequencies}", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = load_frequency_data(args.frequencies)

    if not isinstance(data, dict) or "alleles" not in data:
        print("ERROR: JSON must contain a top-level 'alleles' array", file=sys.stderr)
        return 1
    if not isinstance(data["alleles"], list) or not data["alleles"]:
        print("ERROR: 'alleles' is empty or not an array", file=sys.stderr)
        return 1

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    print(f"Generating figures for {len(data['alleles'])} alleles "
          f"across {data['total_populations']} populations...")

    plot_frequency_heatmap(data, args.output_dir)
    print("  figure1_frequency_heatmap.{png,pdf,svg}")

    plot_gradient_bars(data, args.output_dir)
    print("  figure2_gradient_bars.{png,pdf,svg}")

    plot_population_detail(data, args.output_dir)
    print("  figure3_high_gradient_detail.{png,pdf,svg}")

    print(f"\nAll figures saved to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
