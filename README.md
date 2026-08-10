# lipid_notation

`lipid_notation` converts [Shorthand2020](https://doi.org/10.1194/jlr.S120001025)
lipid names to SMILES/CXSMILES and back. It preserves structural uncertainty
instead of choosing positions or chain assignments that were not measured.

## Quick start

```rust
use lipid_notation::{canonicalize_cxsmiles, name2smiles, smiles2name};

// A fully determined structure needs plain SMILES only.
let determined = name2smiles("FA 18:1(9Z)").unwrap();
assert_eq!(determined, r"OC(=O)CCCCCCC/C=C\CCCCCCCC");

// Unknown double-bond positions are retained in CXSMILES Sg regions.
let ambiguous = name2smiles("FA 18:1").unwrap();
assert!(ambiguous.contains("Sg:"));

// Reverse conversion accepts equivalent atom and branch orders.
let canonical = canonicalize_cxsmiles(&ambiguous).unwrap();
assert_eq!(smiles2name(&canonical).as_deref(), Some("FA 18:1"));
```

Add the crate as a dependency:

```toml
[dependencies]
lipid_notation = { git = "https://github.com/zamboni-lab/lipidoracle-smiles" }
```

The minimum supported Rust version is 1.85.

## Public API

| Function | Purpose |
|---|---|
| `name2smiles(name)` | Convert a lipid name to SMILES or CXSMILES |
| `smiles2name(smiles)` | Canonicalize and convert a supported structure to a lipid name |
| `canonicalize_cxsmiles(smiles)` | Canonicalize SMILES and remap every atom-indexed CX field |
| `name2structure(name)` | Return depiction SMILES plus per-chain atom indexes |
| `expand_cxsmiles_for_depiction(smiles)` | Expand variable-length regions into one drawable representative |
| `class_needs_multi_chain(class)` | Report whether sum composition is structurally ambiguous for a class |

All conversion functions return `None` when a structure cannot be represented
without inventing information. `expand_cxsmiles_for_depiction` returns the input
unchanged when no expansion is necessary.

## How ambiguity is encoded

Plain SMILES is emitted when the name determines one structure. Otherwise the
output uses a small CXSMILES extension model:

| Encoding | Meaning |
|---|---|
| bare `C=C` | Double-bond geometry is unspecified |
| `Sg:` plus `constrain(...)` | A double bond is present in an unlocalized chain region |
| `m:` | A modification is attached to one of several candidate atoms |
| `$snN$` labels plus `swappable(...)` | Explicit chains have unresolved sn assignment |

The `constrain(...)` and `swappable(...)` tokens follow the CXSMILES block in
the SMILES title field. They preserve lipid-specific semantics but are not
standard CXSMILES fields; a toolkit may discard them when serializing again.

`canonicalize_cxsmiles` canonicalizes the molecular graph and updates every
atom reference in `Sg:`, `m:`, atom labels, and atom properties. Lipid trailer
tokens are retained and normalized with their referenced labels. Calling the
function twice returns the same string.

`smiles2name` always canonicalizes first. The reverse converter therefore does
not depend on the exact atom or branch order produced by `name2smiles`. Before
returning a name, it regenerates the structure and verifies canonical CXSMILES
equivalence.

## Lipid notation

Explicit multi-chain names distinguish assignment confidence:

- `/` means the sn positions are resolved, for example `PC 16:0/18:1(9Z)`.
- `_` means the chains are known but their sn assignment is unresolved, for
  example `PC 16:0_18:1`.

Species-level sum compositions such as `PC 34:1`, `DG 36:2`, and `TG 54:3`
are rejected because many chain combinations match them. Single-chain names
such as `FA 18:1` remain valid.

Supported classes include:

- fatty-acyl and neutral-lipid forms: `FA`, `AMP-FA`, `NAE`, `CAR`, `CE`,
  `MG`, `DG`, `TG`, and `CL`;
- glycerophospholipids and lyso forms: `PC`, `PE`, `PS`, `PG`, `PI`, `PA`,
  `LPC`, `LPE`, `LPS`, `LPG`, `LPI`, and `LPA`;
- sphingolipids: `Cer`, `CerP`, `SM`, `HexCer`, `IPC`, `Sph`/`SB`, and `S1P`;
- the chain-free sterol form `ST`.

Chains may carry localized or unlocalized double bonds, geometry, supported
Table 1A substituents, epoxides, and cyclic groups. Generic oxygen counts such
as `;O2` are rejected unless the functional groups are specified.

## Depiction and toolkit behavior

CXSMILES support varies between chemistry toolkits. In particular, renderers
that ignore `Sg:` see only the fixed scaffold and produce a truncated chain.
Use `expand_cxsmiles_for_depiction` or `name2structure` before depiction. The
expanded result is one representative layout; it does not turn an unknown
position into a measurement.

The position-variation `m:` field and atom labels are standard CXSMILES data.
The lipid-specific trailer should be preserved separately if a toolkit rewrites
the string.

## Limitations

- Carbohydrate sequences such as `Gal-Glc-Cer` do not determine glycosidic
  linkage positions and are not converted.
- A functional group on C1 of an acyl chain can change the linkage itself and
  is rejected.
- Confidence percentages are display metadata and are stripped before
  structural conversion.
- `smiles2name` recognizes structures within this crate's supported lipid
  templates; it is not a general structure-to-name engine.
- Some equivalent input names normalize to one spelling during reverse
  conversion, such as `;ep(5)` to `;5Ep`.

## Development

```bash
cargo fmt --check
cargo test
cargo clippy --all-targets --all-features -- -D warnings
```

The test corpus is stored in `testdata/`:

- `name2smiles.csv` pins generated and depiction strings;
- `chains.csv` validates per-chain carbon counts independently;
- `doc_examples.csv` pins documentation examples.

After an intentional output-format change, run
`cargo run --example bless_testdata` and review the corpus diff. For external
toolkit validation, install RDKit and run:

```bash
python scripts/validate_cxsmiles.py
python scripts/validate_cxsmiles.py --cdk
```

The `--cdk` option sends structures to the public CDK Depict service and
therefore requires network access.

## License

MIT
