"""Schema validation and cross-reference checks for the allele catalog.

Validates that african_alleles.json conforms to the expected schema,
checks internal consistency, and flags potential data quality issues.

Usage:
    python validate_catalog.py --catalog /path/to/african_alleles.json
"""

import argparse
import json
import math
import sys
from pathlib import Path

VALID_FUNCTIONS = {"normal_function", "decreased_function", "no_function"}
VALID_EVIDENCE_LEVELS = {"L1", "L2"}
VALID_POPULATIONS = {"YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW"}
VALID_ACTIVITY_SCORES = {0.0, 0.5, 1.0}

REQUIRED_META_FIELDS = {"version", "curation_date", "curator"}
REQUIRED_ENTRY_FIELDS = {
    "gene", "allele_name", "defining_variants", "function",
    "activity_score", "evidence_level", "populations",
    "frequency_range", "source_pmids", "in_pharmvar", "notes",
}
REQUIRED_VARIANT_FIELDS = {"rsid", "position", "ref", "alt"}


def validate_meta(meta: dict, errors: list[str], warnings: list[str]) -> None:
    """Check _meta block for required fields."""
    if not isinstance(meta, dict):
        errors.append("'_meta' must be a JSON object, got " + type(meta).__name__)
        return

    for field in REQUIRED_META_FIELDS:
        if field not in meta:
            errors.append(f"_meta missing required field: {field}")

    version = meta.get("version", "")
    parts = version.split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        warnings.append(f"_meta.version '{version}' is not semver (expected X.Y.Z)")


def validate_entry(idx: int, entry: dict, errors: list[str], warnings: list[str]) -> None:
    """Validate a single catalog entry."""
    if not isinstance(entry, dict):
        errors.append(f"entries[{idx}]: must be a JSON object, got {type(entry).__name__}")
        return

    prefix = f"entries[{idx}] ({entry.get('gene', '?')}/{entry.get('allele_name', '?')})"

    # Required fields
    for field in REQUIRED_ENTRY_FIELDS:
        if field not in entry:
            errors.append(f"{prefix}: missing required field '{field}'")

    # Gene name format
    gene = entry.get("gene", "")
    if gene and not gene[0].isupper():
        warnings.append(f"{prefix}: gene name '{gene}' should be uppercase")

    # Allele name format (should start with *)
    allele = entry.get("allele_name", "")
    if allele and not allele.startswith("*"):
        errors.append(f"{prefix}: allele_name '{allele}' must start with '*'")

    # Function validation
    func = entry.get("function", "")
    if func and func not in VALID_FUNCTIONS:
        errors.append(
            f"{prefix}: function '{func}' not in {sorted(VALID_FUNCTIONS)}"
        )

    # Activity score validation
    score = entry.get("activity_score")
    if score is not None:
        if isinstance(score, bool) or not isinstance(score, (int, float)):
            errors.append(
                f"{prefix}: activity_score must be a number, got {type(score).__name__}"
            )
        elif score not in VALID_ACTIVITY_SCORES:
            warnings.append(
                f"{prefix}: activity_score {score} is non-standard "
                f"(expected one of {sorted(VALID_ACTIVITY_SCORES)})"
            )

    # Activity score consistency with function
    if func == "no_function" and score != 0.0:
        errors.append(f"{prefix}: no_function allele should have activity_score 0.0, got {score}")
    if func == "normal_function" and score != 1.0:
        errors.append(f"{prefix}: normal_function allele should have activity_score 1.0, got {score}")

    # Evidence level
    level = entry.get("evidence_level", "")
    if level and level not in VALID_EVIDENCE_LEVELS:
        errors.append(f"{prefix}: evidence_level '{level}' not in {sorted(VALID_EVIDENCE_LEVELS)}")

    # Populations
    pops = entry.get("populations", [])
    if not pops:
        warnings.append(f"{prefix}: no populations listed")
    for pop in pops:
        if pop not in VALID_POPULATIONS:
            warnings.append(f"{prefix}: population '{pop}' not in expected set")

    # Frequency range format (expected: "X.XX-Y.YY")
    freq_range = entry.get("frequency_range", "")
    if freq_range:
        parts = freq_range.split("-")
        if len(parts) != 2:
            errors.append(f"{prefix}: frequency_range '{freq_range}' must be 'low-high'")
        else:
            try:
                low, high = float(parts[0]), float(parts[1])
                if not math.isfinite(low) or not math.isfinite(high):
                    errors.append(f"{prefix}: frequency_range contains non-finite value")
                else:
                    if low < 0.0:
                        errors.append(f"{prefix}: frequency_range low ({low}) is negative")
                    if high > 1.0:
                        errors.append(f"{prefix}: frequency_range high ({high}) exceeds 1.0")
                    if low > high:
                        errors.append(f"{prefix}: frequency_range low ({low}) > high ({high})")
            except ValueError:
                errors.append(f"{prefix}: frequency_range '{freq_range}' contains non-numeric values")

    # Source PMIDs
    pmids = entry.get("source_pmids", [])
    if not pmids:
        warnings.append(f"{prefix}: no source_pmids — traceability gap")
    for pmid in pmids:
        if not str(pmid).isdigit():
            errors.append(f"{prefix}: PMID '{pmid}' is not numeric")

    # Defining variants
    variants = entry.get("defining_variants", [])
    # CYP3A5*1 is the reference allele, so empty defining_variants is acceptable
    if not variants and allele != "*1":
        warnings.append(f"{prefix}: no defining_variants (only valid for reference alleles)")

    seen_rsids = set()
    for vi, var in enumerate(variants):
        var_prefix = f"{prefix}.defining_variants[{vi}]"
        if not isinstance(var, dict):
            errors.append(f"{var_prefix}: must be a JSON object, got {type(var).__name__}")
            continue

        for field in REQUIRED_VARIANT_FIELDS:
            if field not in var:
                errors.append(f"{var_prefix}: missing '{field}'")

        rsid = var.get("rsid", "")
        if rsid and not rsid.startswith("rs"):
            errors.append(f"{var_prefix}: rsid '{rsid}' must start with 'rs'")
        if rsid in seen_rsids:
            errors.append(f"{var_prefix}: duplicate rsid '{rsid}'")
        seen_rsids.add(rsid)

        pos = var.get("position")
        if pos is not None and (not isinstance(pos, int) or pos <= 0):
            errors.append(f"{var_prefix}: position must be a positive integer, got {pos}")

        ref = var.get("ref", "")
        alt = var.get("alt", "")
        if ref and not all(c in "ACGTN" for c in ref.upper()):
            errors.append(f"{var_prefix}: ref '{ref}' contains invalid bases")
        if alt and not all(c in "ACGTN" for c in alt.upper()):
            errors.append(f"{var_prefix}: alt '{alt}' contains invalid bases")
        if ref and alt and ref == alt:
            errors.append(f"{var_prefix}: ref == alt ('{ref}')")


def validate_cross_references(entries: list[dict], errors: list[str], warnings: list[str]) -> None:
    """Check for duplicates and cross-entry consistency."""
    seen = set()
    for entry in entries:
        key = (entry.get("gene", ""), entry.get("allele_name", ""))
        if key in seen:
            errors.append(f"Duplicate entry: {key[0]} {key[1]}")
        seen.add(key)

    # Check that rsIDs are not shared across different genes
    # (same rsID in different alleles of the same gene is fine)
    rsid_to_genes: dict[str, set[str]] = {}
    for entry in entries:
        gene = entry.get("gene", "")
        for var in entry.get("defining_variants", []):
            rsid = var.get("rsid", "")
            if rsid:
                rsid_to_genes.setdefault(rsid, set()).add(gene)

    for rsid, genes in rsid_to_genes.items():
        if len(genes) > 1:
            warnings.append(
                f"rsID {rsid} appears in multiple genes: {sorted(genes)} "
                f"— verify this is intentional (shared variant)"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate allele catalog JSON")
    parser.add_argument(
        "--catalog", type=Path,
        default=Path(__file__).resolve().parent.parent.parent.parent
        / "data" / "alleles" / "african_alleles.json",
    )
    args = parser.parse_args()

    if not args.catalog.exists():
        print(f"ERROR: catalog not found: {args.catalog}", file=sys.stderr)
        return 1

    with open(args.catalog, encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"ERROR: invalid JSON: {e}", file=sys.stderr)
            return 1

    errors: list[str] = []
    warnings: list[str] = []

    # Structure checks
    if "_meta" not in data:
        errors.append("Missing top-level '_meta' object")
    else:
        validate_meta(data["_meta"], errors, warnings)

    if "entries" not in data:
        errors.append("Missing top-level 'entries' array")
    elif not isinstance(data["entries"], list):
        errors.append("'entries' must be an array")
    else:
        for idx, entry in enumerate(data["entries"]):
            validate_entry(idx, entry, errors, warnings)
        validate_cross_references(data["entries"], errors, warnings)

    # Report
    n_entries = len(data.get("entries", []))
    print(f"Catalog: {args.catalog}")
    print(f"Version: {data.get('_meta', {}).get('version', 'unknown')}")
    print(f"Entries: {n_entries}")
    print()

    if errors:
        print(f"ERRORS ({len(errors)}):")
        for e in errors:
            print(f"  {e}")
        print()

    if warnings:
        print(f"WARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  {w}")
        print()

    if not errors and not warnings:
        print("PASSED — no issues found")
    elif not errors:
        print(f"PASSED with {len(warnings)} warning(s)")
    else:
        print(f"FAILED — {len(errors)} error(s), {len(warnings)} warning(s)")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
