"""Top-level pipeline for all afripharmagen-catalog analyses.

Runs each analysis module in sequence with the correct paths resolved
relative to this file's location. Can be invoked from any working directory.

Usage:
    uv run python run_analyses.py
    uv run python run_analyses.py --pharmcat-dir /path/to/pharmcat_results/reports/
    uv run python run_analyses.py --skip-concordance --skip-plots

Environment:
    Run `uv sync` to install dependencies from pyproject.toml into a local .venv.
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Standard data paths
CATALOG = REPO_ROOT / "data" / "alleles" / "african_alleles.json"
SAMPLES = REPO_ROOT / "data" / "samples" / "1000g_african_samples.tsv"
BENCHMARK = REPO_ROOT / "results" / "afripharmagen_benchmark.jsonl"

# Analysis scripts
FREQ_SCRIPT = REPO_ROOT / "analyses" / "population_frequencies" / "scripts" / "estimate_frequencies.py"
HWE_SCRIPT = REPO_ROOT / "analyses" / "population_frequencies" / "scripts" / "hwe_test.py"
PLOT_SCRIPT = REPO_ROOT / "analyses" / "population_frequencies" / "scripts" / "plot_gradients.py"
CONCORDANCE_SCRIPT = REPO_ROOT / "analyses" / "concordance_benchmark" / "scripts" / "compare_pharmcat.py"
VALIDATE_SCRIPT = REPO_ROOT / "analyses" / "allele_catalog" / "scripts" / "validate_catalog.py"

# Output directories
FREQ_OUTPUT = REPO_ROOT / "analyses" / "population_frequencies" / "results"
CONCORDANCE_OUTPUT = REPO_ROOT / "analyses" / "concordance_benchmark" / "results"


def run_step(name: str, cmd: list[str]) -> bool:
    """Run a subprocess step, printing status. Returns True on success.

    All command arguments are constructed internally from resolved Path objects
    and sys.executable — no external/user-controlled strings are interpolated.
    """
    print(f"\n{'─' * 60}")
    print(f"  {name}")
    print(f"{'─' * 60}\n")
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))  # noqa: S603
    if result.returncode != 0:
        print(f"\n  FAILED (exit {result.returncode})", file=sys.stderr)
        return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run all afripharmagen-catalog analyses end-to-end."
    )
    parser.add_argument(
        "--pharmcat-dir", type=Path, default=None,
        help="Path to PharmCAT JSON report directory. Skip concordance if omitted.",
    )
    parser.add_argument(
        "--skip-concordance", action="store_true",
        help="Skip the PharmCAT concordance benchmark.",
    )
    parser.add_argument(
        "--skip-plots", action="store_true",
        help="Skip gradient plot generation (useful in headless environments).",
    )
    args = parser.parse_args()

    python = sys.executable
    failed = []

    # Preflight: check data files exist
    for path, label in [
        (CATALOG, "allele catalog"),
        (SAMPLES, "sample metadata"),
        (BENCHMARK, "benchmark JSONL"),
    ]:
        if not path.exists():
            print(f"ERROR: {label} not found at {path}", file=sys.stderr)
            return 1

    # Step 1: Validate catalog schema
    if VALIDATE_SCRIPT.exists():
        ok = run_step(
            "Catalog Validation",
            [python, str(VALIDATE_SCRIPT), "--catalog", str(CATALOG)],
        )
        if not ok:
            failed.append("catalog validation")
    else:
        print(f"\n  SKIP: {VALIDATE_SCRIPT.name} not found")

    # Step 2: Population frequency estimation
    ok = run_step(
        "Population Frequency Estimation",
        [
            python, str(FREQ_SCRIPT),
            "--benchmark", str(BENCHMARK),
            "--catalog", str(CATALOG),
            "--samples", str(SAMPLES),
            "--output-dir", str(FREQ_OUTPUT),
        ],
    )
    if not ok:
        failed.append("frequency estimation")

    # Step 3: Hardy-Weinberg equilibrium test
    ok = run_step("Hardy-Weinberg Equilibrium Test", [python, str(HWE_SCRIPT)])
    if not ok:
        failed.append("HWE test")

    # Step 4: Frequency gradient plots
    if not args.skip_plots and PLOT_SCRIPT.exists():
        ok = run_step(
            "Frequency Gradient Plots",
            [
                python, str(PLOT_SCRIPT),
                "--frequencies", str(FREQ_OUTPUT / "population_frequencies.json"),
                "--output-dir", str(FREQ_OUTPUT),
            ],
        )
        if not ok:
            failed.append("gradient plots")
    elif args.skip_plots:
        print("\n  SKIP: plots (--skip-plots)")
    else:
        print(f"\n  SKIP: {PLOT_SCRIPT.name} not found")

    # Step 5: PharmCAT concordance benchmark
    if args.skip_concordance or args.pharmcat_dir is None:
        print("\n  SKIP: concordance benchmark (no --pharmcat-dir provided)")
    else:
        # Preflight: verify directory exists and contains reports
        if not args.pharmcat_dir.is_dir():
            print(
                f"ERROR: --pharmcat-dir is not a directory: {args.pharmcat_dir}",
                file=sys.stderr,
            )
            failed.append("concordance benchmark")
        elif not list(args.pharmcat_dir.glob("*.report.json")):
            print(
                f"ERROR: no *.report.json files found in {args.pharmcat_dir}",
                file=sys.stderr,
            )
            failed.append("concordance benchmark")
        else:
            ok = run_step(
                "PharmCAT Concordance Benchmark",
                [
                    python, str(CONCORDANCE_SCRIPT),
                    "--pharmcat-dir", str(args.pharmcat_dir),
                    "--afripharmagen", str(BENCHMARK),
                    "--samples", str(SAMPLES),
                    "--catalog", str(CATALOG),
                    "--output-dir", str(CONCORDANCE_OUTPUT),
                ],
            )
            if not ok:
                failed.append("concordance benchmark")

    # Summary
    print(f"\n{'═' * 60}")
    if failed:
        print(f"  PIPELINE FINISHED WITH FAILURES: {', '.join(failed)}")
        print(f"{'═' * 60}")
        return 1
    print("  PIPELINE COMPLETE — all steps succeeded")
    print(f"{'═' * 60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
