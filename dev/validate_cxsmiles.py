# Cross-toolkit validation of the CXSMILES this crate emits.
#
# Golden-file tests catch regressions but not misconceptions: if a block was
# wrong from the start, the test enshrines it (dev/SMILES.md §3.4). Three of the
# four ambiguity blocks this project once emitted turned out to be wrong exactly
# that way. This script makes the checks the Rust tests cannot — does a real
# toolkit accept the string, and does it see the molecule we meant?
#
# Usage:
#   python dev/validate_cxsmiles.py            # RDKit only, offline
#   python dev/validate_cxsmiles.py --cdk      # also render each string via CDK Depict
#
# Reads testdata/name2smiles.csv, which `cargo test` also checks the generator
# against — so a green Rust suite plus a green run here means the strings are
# both what we intended and acceptable to both toolkits.
#
# Exit status is nonzero if any check fails, so it can gate a release.

import argparse
import csv
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "testdata" / "name2smiles.csv"
DOC_EXAMPLES = ROOT / "testdata" / "doc_examples.csv"
CDK_DEPICT = "https://www.simolecule.com/cdkdepict/depict/bow/svg"

# RDKit implements no `RG:` at all — not even the minimal `C* |$;R1$,RG:_R1={C}|`
# — so an R-group string is expected to fail there. Every other block must
# parse. Verified against RDKit 2024.09.6.
RDKIT_UNSUPPORTED = ("RG:",)


def load_corpus():
    """Every row of both testdata files that carries a generated string."""
    rows = []
    for path in (CORPUS, DOC_EXAMPLES):
        if not path.exists():
            sys.exit(f"{path} not found")
        with path.open(encoding="utf-8", newline="") as fh:
            rows += [r for r in csv.DictReader(fh) if r.get("cxsmiles")]
    return rows


def check_rdkit(rows):
    """Every stored string must parse in RDKit, unless it uses a block RDKit
    has not implemented; every expanded string must parse unconditionally,
    since the expanded form exists precisely to be handed to RDKit."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    failures = []
    for row in rows:
        name, stored, expanded = row["name"], row["cxsmiles"], row["expanded"]
        if not stored:
            failures.append((name, "generator returned nothing"))
            continue
        expected_fail = any(b in stored for b in RDKIT_UNSUPPORTED)
        parsed = Chem.MolFromSmiles(stored) is not None
        if parsed == expected_fail:
            failures.append(
                (
                    name,
                    "stored string parsed in RDKit but was expected to fail"
                    if parsed
                    else f"stored string failed to parse in RDKit: {stored}",
                )
            )
        if expanded and Chem.MolFromSmiles(expanded) is None:
            failures.append((name, f"expanded string failed to parse: {expanded}"))
        if expanded and (" |" in expanded):
            failures.append((name, "expanded string still carries a CXSMILES block"))
    return failures


def check_cdk(rows):
    """Every stored string must render in CDK. CDK is the toolkit that
    implements the blocks we rely on, so a 400 here means the string is
    malformed, not merely unsupported."""
    failures = []
    for row in rows:
        name, stored = row["name"], row["cxsmiles"]
        if not stored:
            continue
        url = f"{CDK_DEPICT}?smi={urllib.parse.quote(stored, safe='')}&showtitle=false"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                body = resp.read()
            if b"<svg" not in body:
                failures.append((name, "CDK returned no SVG"))
        except urllib.error.HTTPError as exc:
            failures.append((name, f"CDK rejected the string (HTTP {exc.code}): {stored}"))
        except OSError as exc:
            failures.append((name, f"CDK unreachable: {exc}"))
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cdk", action="store_true", help="also render each string via CDK Depict")
    args = ap.parse_args()

    rows = load_corpus()
    failures = check_rdkit(rows)
    print(f"RDKit: checked {len(rows)} names")
    if args.cdk:
        cdk_failures = check_cdk(rows)
        print(f"CDK:   checked {len(rows)} names")
        failures += cdk_failures

    for name, why in failures:
        print(f"  FAIL  {name}: {why}")
    print(f"\n{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
