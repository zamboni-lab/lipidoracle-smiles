# lipid_notation

Convert [Shorthand2020](https://doi.org/10.1194/jlr.S120001025) lipid names to
SMILES/CXSMILES, encoding **what was and wasn't determined** instead of guessing
a fully specified structure.

Extracted from [LipidOracle](https://github.com/zamboni-lab/lipidoracle) so it
can be used, reviewed and cited on its own. No dependencies.

```rust
use lipid_notation::name2smiles;

// Fully determined -> plain, ordinary SMILES.
name2smiles("FA 18:1(9Z)");   // Some(r"OC(=O)CCCCCCC/C=C\CCCCCCCC")

// Double bond present, position unknown -> a variable-length run plus a size
// constraint, never a guessed position.
name2smiles("FA 18:1");       // Some("OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)")

// Hydroxyl present, position unknown -> a position-variation bond over every
// candidate carbon.
name2smiles("FA 18:0;OH");    // Some("OC(=O)CCCCCCCCCCCCCCCCC.*O |m:20:3.4.…19|")

// sn-position unknown -> linking atoms labelled, assignment left permutable.
name2smiles("DG 16:0_18:1(9)");
// Some("C(CO)(OC(=O)CCCCCCCC=CCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC //        |$;;;sn2;;;;;;;;;;;;;;;;;;;;;sn1$| swappable(sn1,sn2)")

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
| `$snN$` + `swappable(...)` | "either chain could be at either position" | chains joined with `_` |

`dev/SMILES.md` §3 records what was wrong, why it survived a full test suite,
and what replaced it. Corrections came from review by John Mayfield.

## API

| function | status |
|---|---|
| `name2smiles(name) -> Option<String>` | the stored CXSMILES |
| `name2structure(name) -> Option<LipidStructure>` | plain SMILES + per-chain atom indices, for highlighting |
| `expand_cxsmiles_for_depiction(smi) -> String` | stored → drawable by any tool |
| `class_needs_multi_chain(class) -> bool` | whether a class rejects sum-composition shorthand |
| `smiles2name(smi) -> Option<String>` | the stored string back to a name |

`smiles2name` inverts `name2smiles` — it reads this crate's own output, not
arbitrary third-party SMILES. Every answer is proved before it is returned: the
name is fed back through `name2smiles` and must regenerate the input exactly, so
the function can lose coverage but cannot return a name meaning something other
than the structure given. That is what makes the corpus round-trip test a
validation of the forward direction rather than a formatting check — see
`dev/SMILES.md` §4.6.

## Toolkit support

Stored strings are rigorous, not universally drawable. **Verified against RDKit
2024.09.6 and CDK Depict:**

| block | CDK | RDKit |
|---|---|---|
| plain `C=C` | ✅ | ✅ |
| `Sg:` + size constraint | ✅ labelled repeat brackets | ❌ **silently ignored** |
| `m:` (`*O` dummy-stub form) | ✅ collapses to one position | ✅ highlights all candidates |
| `$snN$` labels + `swappable(...)` | ✅ atom labels | ✅ labels survive canonicalization |

The RDKit `Sg:` row is the dangerous one — no warning, no error, just a molecule
missing most of its chain. **Call `expand_cxsmiles_for_depiction` before handing
a stored string to RDKit.**

The last row used to read `RG:` / ❌ hard parse failure. Replacing it with atom
labels plus a trailing token means every string this crate emits now parses in
both toolkits. The trailing token itself is *not* CXSMILES — it rides in the
SMILES title field and is dropped when a toolkit re-serializes the molecule, so
it must never be the only place a caveat lives. `dev/extension.md` §5.1.

## Testing

```bash
cargo test                              # unit, golden and property tests
```

A cross-toolkit validator (`validate_cxsmiles.py`, RDKit and CDK Depict) lives
with the design notes rather than in this repository — see *Documentation*.

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
| `doc_examples.csv` | strings quoted in the design notes, pinned so the prose cannot drift |

## Known limitations

Documented at length in `dev/SMILES.md` §4:

- **Table 1C carbohydrates are not handled.** `Gal-Glc-Cer` names the sugar
  sequence but never the glycosidic linkage positions, so there is no single
  honest structure to emit.
- **A group on C1 of an acyl chain is refused.** `FA 18:0;1OMe` is the methyl
  ester — a different linkage, not a substituted carbon, and the chain builder
  cannot change the headgroup it is attached to.
- **Confidence percentages have nowhere to live.** No CXSMILES construct
  expresses *weighted* alternatives, so a 51% call and a 100% call are written
  identically.
- **Variable names run out at ten**, so a cardiolipin with four unlocalized
  chains emits a constraint set that is inconsistent read as algebra.
- **Species-level shorthand is rejected** rather than served, which is the most
  common case in published data.

## Documentation

The `dev/` design notes are kept outside this repository. Sections of them are
cited throughout this README and in the rustdoc, and the citations are stable:

- `SMILES.md` — the technical write-up: every block, every trade-off, every
  open problem, with depictions from both toolkits.
- `NOMENCLATURE.md` — the supported name grammar, class by class.
- `extension.md` — why the trailing token language replaced CXSMILES `RG:`.
- `validate_cxsmiles.py` — the cross-toolkit validator.

## License

MIT
