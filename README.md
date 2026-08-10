# shorthand2smiles

Convert [Shorthand2020](https://doi.org/10.1194/jlr.S120001025) lipid names to
SMILES/CXSMILES, encoding **what was and wasn't determined** instead of guessing
a fully specified structure.

Extracted from [LipidOracle](https://github.com/zamboni-lab/lipidoracle) so it
can be used, reviewed and cited on its own. No dependencies.

```rust
use shorthand2smiles::name2smiles;

// Fully determined -> plain, ordinary SMILES.
name2smiles("FA 18:1(9Z)");   // Some(r"OC(=O)CCCCCCC/C=C\CCCCCCCC")

// Double bond present, position unknown -> a variable-length run plus a size
// constraint, never a guessed position.
name2smiles("FA 18:1");       // Some("OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| a+b=15")

// Hydroxyl present, position unknown -> a position-variation bond over every
// candidate carbon.
name2smiles("FA 18:0;OH");    // Some("OC(=O)CCCCCCCCCCCCCCCCC.*O |m:20:3.4.…19|")

// sn-position unknown -> the chains become interchangeable R-group alternatives.
name2smiles("DG 16:0_18:1(9)");
// Some("C(CO)(O*)CO* |$;;;;R1;;;R1$,RG:_R1={C(=O)CCCCCCCCCCCCCCC},{C(=O)CCCCCCCC=CCCCCCCCC}|")

// Sum composition only -> nothing, because any single answer would be invented.
name2smiles("PC 34:1");       // None
```

The presence of a `|...|` tail is itself the signal: **pipes mean something was
undetermined.**

## Why this exists

Lipidomics reports structure at whatever level the evidence supports. `PC 34:1`
concedes everything but the sum composition; `PC 16:0_18:1` names the chains but
not their positions; `PC 16:0/18:1(9Z)` concedes nothing. SMILES only speaks the
last of those, so the usual options are to invent the missing detail or emit
nothing at all. This crate takes the third option.

| block | says | emitted when |
|---|---|---|
| *(none)* | geometry undetermined | always — a bare `C=C` already means "cis or trans, not determined" |
| `Sg:` | "the double bond is somewhere in this stretch" | a chain declares more double bonds than it localizes |
| `m:` | "this group attaches somewhere on this chain" | a modification is written with no position (`;OH`) |
| `RG:` | "either chain could be at either position" | chains joined with `_`, and none needs `Sg:`/`m:` |

Two CXSMILES constructs are deliberately **not** used, both of which earlier
revisions of this code got wrong:

- **`ctu:`** is a *query* feature for matching either configuration when
  searching. A plain `C=C` already says the same thing about a structure, so the
  block added nothing and made every partially-determined structure a query.
- **`f:`** groups components into one entity (a salt, a hydrate). It expresses
  *and*, never *or*, so it cannot say "this chain **or** that chain sits at
  sn-1".

`dev/SMILES.md` §3 records what was wrong, why it survived a full test suite,
and what replaced it. Corrections came from review by John Mayfield.

## API

| function | status |
|---|---|
| `name2smiles(name) -> Option<String>` | the stored CXSMILES |
| `name2structure(name) -> Option<LipidStructure>` | plain SMILES + per-chain atom indices, for highlighting |
| `expand_cxsmiles_for_depiction(smi) -> String` | stored → drawable by any tool |
| `class_needs_multi_chain(class) -> bool` | whether a class rejects sum-composition shorthand |
| `smiles2name(smi) -> Option<String>` | **not implemented yet** — always `None` |

`smiles2name` has a fixed signature and a written, `#[ignore]`d round-trip test,
so implementing it should disturb no caller. It is worth building both as
validation (the forward mapping has no parser-level check that a string means
what the name meant) and to allow ingesting structures from other tools.

## Toolkit support

Stored strings are rigorous, not universally drawable. **Verified against RDKit
2024.09.6 and CDK Depict:**

| block | CDK | RDKit |
|---|---|---|
| plain `C=C` | ✅ | ✅ |
| `Sg:` + size constraint | ✅ labelled repeat brackets | ❌ **silently ignored** |
| `m:` (`*O` dummy-stub form) | ✅ collapses to one position | ✅ highlights all candidates |
| `RG:` + `$...$` labels | ✅ Markush scheme | ❌ hard parse failure |

The RDKit `Sg:` row is the dangerous one — no warning, no error, just a molecule
missing most of its chain. **Call `expand_cxsmiles_for_depiction` before handing
a stored string to RDKit.**

## Testing

```bash
cargo test                              # unit, golden and property tests
python dev/validate_cxsmiles.py         # RDKit parse checks, offline
python dev/validate_cxsmiles.py --cdk   # also render every string via CDK Depict
```

Two kinds of test, and the distinction is the point:

- **Golden** — `testdata/name2smiles.csv` records the exact strings emitted.
  Catches regressions. Cannot catch a misconception: if an expected string was
  wrong from the start, comparing against it passes forever, which is exactly
  how three wrong blocks survived for as long as they did.
- **Property** — `testdata/chains.csv` holds per-chain carbon counts read off
  the *names* by hand, and the integration tests assert structural invariants
  (every `m:` target is a `*`, every `R1` label lands on a wildcard, no block is
  a query feature). These are statements about what the output must satisfy, not
  recordings of what it currently does.

After a deliberate encoding change, `cargo run --example bless_testdata` rewrites
the golden columns; then read `git diff testdata/` carefully, because that diff
is the only review a golden file ever gets. `chains.csv` is never regenerated.

### testdata

| file | contents |
|---|---|
| `name2smiles.csv` | `name, cxsmiles, expanded` — 54 names covering every shape of ambiguity the generator emits |
| `chains.csv` | `name, sn, carbons` — hand-written per-chain carbon counts |
| `doc_examples.csv` | strings quoted in `dev/`, including one documenting a known bug |

## Known limitations

Documented at length in `dev/SMILES.md` §4:

- **Rings are silently dropped.** `FA 18:0;ep(5)` returns plain stearic acid
  rather than an epoxystearate — a different molecule, reported as fact. The one
  outright bug rather than a trade-off.
- **Unresolved sn-position is lost whenever a chain needs `Sg:`.** CDK rejects a
  nested block inside an `RG:` definition, so `PC 16:0_18:1` can express the
  double-bond ambiguity or the sn ambiguity but not both. It keeps the former.
- **Confidence percentages have nowhere to live.** No CXSMILES construct
  expresses *weighted* alternatives, so a 51% call and a 100% call are written
  identically.
- **Variable names run out at ten**, so a cardiolipin with four unlocalized
  chains emits a constraint set that is inconsistent read as algebra.
- **Species-level shorthand is rejected** rather than served, which is the most
  common case in published data.

## Documentation

- `dev/SMILES.md` — the technical write-up: every block, every trade-off, every
  open problem, with depictions from both toolkits.
- `dev/NOMENCLATURE.md` — the supported name grammar, class by class.
- `dev/validate_cxsmiles.py` — the cross-toolkit validator.

## License

MIT
