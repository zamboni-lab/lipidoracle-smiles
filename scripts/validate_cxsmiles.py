# Cross-toolkit validation of the CXSMILES corpus.
#
# Usage:
#   python scripts/validate_cxsmiles.py        # RDKit only, offline
#   python scripts/validate_cxsmiles.py --cdk  # also render via CDK Depict
#
# Exit status is nonzero when a stored or expanded structure is rejected.

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
    """Every stored string must parse in RDKit, and every expanded string must
    too, since the expanded form exists precisely to be handed to RDKit.

    The trailer must survive as the title field, and its `$snN$` anchors must
    survive canonicalization so `swappable(...)` can still name them."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    failures = []
    for row in rows:
        name, stored, expanded = row["name"], row["cxsmiles"], row["expanded"]
        if not stored:
            failures.append((name, "generator returned nothing"))
            continue
        parsed = Chem.MolFromSmiles(stored) is not None
        if not parsed:
            failures.append((name, f"stored string failed to parse in RDKit: {stored}"))
        if expanded and Chem.MolFromSmiles(expanded) is None:
            failures.append((name, f"expanded string failed to parse: {expanded}"))
        if expanded and (" |" in expanded):
            failures.append((name, "expanded string still carries a CXSMILES block"))

        failures += check_trailer(Chem, name, stored)
    return failures


def check_trailer(Chem, name, stored):
    """The trailer is our token list, and it rides in the SMILES title field.

    Two properties are checked beyond successful parsing:

    1. The whole trailer arrives intact — parentheses, commas and semicolons
       included — or a token silently loses its arguments.
    2. Any `snN` atom label survives canonical rewriting *and stays on its own
       atom*, allowing `swappable()` to refer to stable labels instead of atom
       indexes.
    """
    if " |" not in stored:
        return []
    trailer = stored.split("|")[-1].strip()
    if not trailer:
        return []

    out = []
    mol = Chem.MolFromSmiles(stored)
    if mol is None:
        return []
    got = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    if got != trailer:
        out.append((name, f"trailer did not survive parsing: {got!r} != {trailer!r}"))

    labels = sorted(
        a.GetProp("atomLabel") for a in mol.GetAtoms() if a.HasProp("atomLabel")
    )
    named = sorted(
        n for t in trailer.split(";")
        if t.startswith("swappable(")
        for n in t[len("swappable("):].rstrip(")").split(",")
    )
    if named and labels != named:
        out.append((name, f"swappable names {named} but the molecule carries {labels}"))

    if labels:
        again = Chem.MolFromSmiles(Chem.MolToCXSmiles(mol))
        survived = sorted(
            a.GetProp("atomLabel") for a in again.GetAtoms() if a.HasProp("atomLabel")
        )
        if survived != labels:
            out.append(
                (name, f"labels lost on canonical rewrite: {labels} -> {survived}")
            )
    return out


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
