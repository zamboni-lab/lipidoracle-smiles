#!/usr/bin/env python3
"""Validate the CXSMILES examples in blog/v2_followup.md.

Backs the post's claim that every string parses in both CDK and RDKit and
that the lipid trailer survives as the SMILES title field.

Usage:
    python blog/validate_blog_examples.py          # RDKit only, offline
    python blog/validate_blog_examples.py --cdk    # also render via CDK Depict

Exit status is nonzero when any example is rejected.

Run inside a Python that has RDKit installed, e.g.
    uv run --python 3.12 --with rdkit python blog/validate_blog_examples.py --cdk
"""
import argparse
import re
import urllib.parse
import urllib.request

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "blog" / "v2_followup.md"


def extract_pairs(md: str):
    """Pull every `NAME -> CXSMILES` pair out of the ```text blocks."""
    pairs = []
    blocks = re.findall(r"```text\n(.*?)```", md, re.S)
    for b in blocks:
        for chunk in re.split(r"\n(?=\S.*\u2192)", b.strip()):
            chunk = chunk.strip()
            if "\u2192" not in chunk:
                continue
            name, smi = chunk.split("\u2192", 1)
            pairs.append((name.strip(), smi.strip()))
    return pairs


def check_rdkit(pairs):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    fail = []
    for name, smi in pairs:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fail.append((name, "RDKit failed to parse"))
            continue
        if " |" in smi:
            trailer = smi.split("|")[-1].strip()
            prop = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
            if trailer and prop != trailer:
                fail.append((name, f"trailer lost in _Name: got {prop!r}"))
    return fail


def check_cdk(pairs):
    fail = []
    for name, smi in pairs:
        enc = urllib.parse.quote(smi, safe="-_.~")
        url = "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=" + enc
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            r = urllib.request.urlopen(req, timeout=90)
            data = r.read()
            ct = r.headers.get("Content-Type", "")
            if not (r.status == 200 and ("image/svg+xml" in ct or b"<svg" in data[:2000])):
                fail.append((name, f"CDK rejected (status {r.status}, {ct})"))
        except Exception as e:  # noqa: BLE001
            fail.append((name, f"CDK error: {e}"))
    return fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdk", action="store_true", help="also validate via CDK Depict")
    args = ap.parse_args()

    md = MD.read_text(encoding="utf-8")
    pairs = extract_pairs(md)
    print(f"blog examples: {len(pairs)}")

    fail = check_rdkit(pairs)
    print(f"RDKit: {len(pairs) - len(fail)}/{len(pairs)} passed")
    if args.cdk:
        fail2 = check_cdk(pairs)
        print(f"CDK  : {len(pairs) - len(fail2)}/{len(pairs)} passed")
        fail = fail + fail2

    for name, msg in fail:
        print(f"  FAIL {name}: {msg}")
    if fail:
        print("RESULT: FAIL")
        raise SystemExit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
