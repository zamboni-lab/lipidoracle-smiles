# lipid_notation

`lipid_notation` converts [Shorthand2020](https://doi.org/10.1194/jlr.S120001025)
lipid names to SMILES/CXSMILES and back. It preserves structural uncertainty
instead of choosing positions or chain assignments that were not measured.

## Quick start

```rust
use lipid_notation::{canonicalize, name2smiles, smiles_for_depiction, smiles2name};

// A fully determined structure needs plain SMILES only.
let determined = name2smiles("FA 18:1(9Z)").unwrap();
assert_eq!(determined, r"OC(=O)CCCCCCC/C=C\CCCCCCCC");

// Unknown double-bond positions are retained in CXSMILES Sg regions.
let ambiguous = name2smiles("FA 18:1").unwrap();
assert!(ambiguous.contains("Sg:"));

// Reverse conversion accepts equivalent atom and branch orders.
let canonical = canonicalize(&ambiguous).unwrap();
assert_eq!(smiles2name(&canonical).as_deref(), Some("FA 18:1"));

// Make the attachment chosen for an unlocalized group deterministic in a drawing.
let drawing = smiles_for_depiction(&ambiguous);
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
| `canonicalize(smiles)` | Canonicalize SMILES and remap every atom-indexed CX field |
| `smiles_for_depiction(smiles)` | Canonicalize/reindex CXSMILES and select a lipid-friendly bond for each `m:` group |
| `name2structure(name)` | Return depiction SMILES plus per-chain atom indexes |
| `smiles_expand(smiles)` | Expand variable-length regions into one drawable representative |
| `class_needs_multi_chain(class)` | Report whether sum composition is structurally ambiguous for a class |

All conversion functions return `None` when a structure cannot be represented
without inventing information. `smiles_expand` returns the input
unchanged when no expansion is necessary.

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
substituents (according to Liebeisch et al, 2020), epoxides, and cyclic groups. Generic oxygen counts such
as `;O2` are rejected unless the functional groups are specified.

## How ambiguity is encoded in CXSMILES

### The problem

Lipidomics routinely reports structures it has only partly determined, and
SMILES has no way to say so: it describes exactly one molecule, with every bond
order and every substituent position committed. Writing a lipid name into SMILES
therefore means inventing the parts the instrument did not measure.

[CXSMILES](https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html)
recovers some of that, but not all of it:

| What the name says | Can SMILES say it? | Can CXSMILES say it? |
|---|---|---|
| geometry undetermined (`18:1(9)`) | yes — a bare `C=C` already means it | — |
| a double bond somewhere in a stretch (`18:1`) | no | partly: `Sg:` marks the stretch but has nowhere to record its length |
| a modification on one of many carbons (`;OH`) | no | yes: `m:` position-variation |
| chains known, sn assignment not (`16:0_18:1`) | no | badly: `RG:` over-generates and cannot coexist with `Sg:` |
| a weighted call (`Δ9 92%`) | no | no — it goes in the trailer instead |
| sum composition (`PC 34:1`) | no | no |

### The design rule

**Everything the standard covers goes inside the pipes, correct and complete.
Everything it does not goes after them, where no toolkit can mistake it for
structure.**

The part before and inside `|...|` is ordinary, conformant CXSMILES. A toolkit
that knows nothing about lipids parses it, renders it, and gets a chemically
valid molecule. Nothing lipid-specific is smuggled into a standard field, and no
standard field is given a private meaning.

### Inside the pipes: standard CXSMILES

| Encoding | Meaning |
|---|---|
| bare `C=C` | double-bond geometry undetermined — SMILES already means this, so nothing is added |
| `Sg:` | a double bond lies somewhere in this stretch |
| `m:` | this group attaches to one of these candidate atoms |
| `$snN$` atom labels | this atom is the sn-*N* attachment point |

### After the pipes: this crate's tokens

The trailer is a `;`-separated list of `name(argument)` tokens. It generalizes
the `a+b=15` size constraint that the CDK cookbook already puts there, into
something that can carry a second kind of statement:

| Token | Meaning |
|---|---|
| `constrain(a+b=15)` | the `Sg:` runs marked `a` and `b` span 15 carbons between them |
| `swappable(sn1,sn2)` | the chains at these labelled positions may be exchanged |
| `dbPos(sn1:9@92)` | the double bond at Δ9 on sn-1 was called with 92% confidence |
| `mPos(OH1:11OH@50,13OH@50)` | the group on the stub labelled `OH1` is at position 11 or 13, evenly split |

`dbPos` and `mPos` carry the bracketed consensus tail that instrument software
puts after a name — `FA 18:2(9,12) [DB sn1: Δ9 100%, Δ12 88%]`. No structure
format can hold a weighted call, so this used to be stripped and thrown away.
It is metadata, and the trailer is where this crate's metadata lives, so it is
carried there instead and `smiles2name` reconstructs the original tail from it.

Within a token, `,` separates positions and `|` separates *mutually exclusive*
candidates for one feature: `dbPos(sn1:5@100|14@50|15@50)` reads as "Δ5 for
certain, plus one more double bond that is either Δ14 or Δ15". Positions and
percentages are written `pos@percent` with no `Δ` and no spaces, because a
`.smi` reader splits the line on whitespace and would truncate the token; the
original spelling comes back on the way out.

These tokens only ever *refine* what the structure already says. An entry for a
position the SMILES commits to records how sure that call was. A set of `|`
alternatives corresponds to a double bond left inside an `Sg:` run, or a group
left on an `m:` stub — narrowing "somewhere in this stretch" to "one of these,
with these odds" is more information, and still not a determination.

**These tokens are not official CXSMILES.** They sit in what a SMILES reader
treats as the *title* field, so a toolkit reads them as the molecule's name and
never as structure — CDK Depict will print them under the picture unless asked
not to, and RDKit stores them in `_Name` and drops them when it writes the
molecule out again. Strip the trailer and what remains is still valid CXSMILES
describing a chemically valid molecule; you lose the lipid semantics, not the
structure.

Two rules make the trailer safe to extend:

- **Tokens name things, never positions.** `constrain` names `Sg:` variables,
  and `swappable`, `dbPos` and `mPos` name `$...$` atom labels — a chain's
  `snN`, or an `m:` stub's own label. Toolkits maintain those labels across
  canonicalization: renumber the molecule and the label travels with its atom.
  An atom *index* in the trailer would silently rot the first time that
  happened, because nothing rewrites the title field.
- **Anything stated in the trailer is anchored in the pipes.** A `swappable`,
  `dbPos` or `mPos` token always comes with the `$...$` labels it names, so the
  `|...|` block is never empty while the trailer has something to say. That keeps the one-character
  check honest: **pipes mean something was undetermined.**

The same shape extends to anything else worth stating — the grammar is open, and
unknown tokens are ignored rather than treated as errors.

### The roadblocks that forced this

Each of these is a thing CXSMILES cannot express, or expresses wrongly, and each
one shaped the encoding above:

- **An `RG:` definition cannot contain a nested `|...|` block.** CDK rejects it.
  So a chain needing an `Sg:` or `m:` block of its own could never be an R-group
  alternative, and `PC 16:0_18:1` had to choose between saying its double bond
  was unlocalized and saying its sn assignment was. It kept the double bond and
  dropped the sn ambiguity silently. This is why chains are now always written
  into the main string and the sn statement moved to the trailer.
- **A Markush R-group states the wrong cardinality.** Under the ordinary reading
  — one label with *N* definitions substitutes each of its sites independently —
  `PC 16:0_18:1` also admitted `PC 16:0/16:0`, and a TG permitted 27 assignments
  where the name allows 6. Whether a given toolkit reads it that way or pairs
  definitions to sites positionally is itself the problem: nothing in the string
  says which is meant. Labelled positions plus `swappable` are unambiguous.
- **RDKit does not implement `RG:` at all** — not even the minimal case. Dropping
  it means every string this crate emits now parses in both toolkits.
- **`Sg:n` has nowhere to record a repeat count.** Hence `constrain(...)`.
- **`m:` candidate lists cannot see inside `Sg:` runs.** A chain with both an
  unlocalized modification and unlocalized double bonds can only offer the atoms
  physically present in the string, so the candidate list understates the
  ambiguity — the safer direction, but still not the truth.
- **`ctu:` and `f:` look applicable and are not.** `ctu:` is a *query* feature
  for matching either geometry, and a bare `C=C` already says as much about a
  structure; `f:` groups components into one entity and expresses *and*, never
  the *or* that unresolved regiochemistry needs. Both were emitted by earlier
  revisions of this crate and both were wrong.
- **No construct expresses weighted alternatives**, so a 92% call and a 100%
  call are written identically in the structure. The percentages are not lost:
  they move to `dbPos`/`mPos` trailer tokens, which name the chain or the `m:`
  stub they refer to and never touch a CXSMILES field.
- **No construct varies a bond order within a structure**, which is why sum
  compositions such as `PC 34:1` are rejected rather than guessed.

### Examples

Fully determined, so plain SMILES with no tail at all:

```
FA 18:1(9Z)  →  OC(=O)CCCCCCC/C=C\CCCCCCCC
```

One double bond, position not determined — a flexible run plus its length,
never a guessed position:

```
FA 18:1  →  OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)
```

A hydroxyl known to be present but not placed, as a disconnected stub whose
`m:` block lists every carbon it could sit on:

```
FA 18:0;OH  →  OC(=O)CCCCCCCCCCCCCCCCC.*O |m:20:3.4.…19|
```

Both chains known, their sn assignment not:

```
PC 16:0_18:1(9)
→ C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCCCCCCC=CCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC |$;;…;sn2;;…;sn1$| swappable(sn1,sn2)
```

And the case that motivated the whole design — an unlocalized double bond *and*
an unresolved sn assignment, which the `RG:` encoding could not state together:

```
PC 16:0_18:1
→ C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC)COC(=O)CCCCCCCCCCCCCCC |$;;;;;;;;;;;;;sn2;;;;;;;;sn1$,Sg:n:16:a:ht,Sg:n:19:b:ht| swappable(sn1,sn2);constrain(a+b=15)
```

Worth noting: each chain gets its own equation, and equations are matched to
`Sg:` markers *positionally*, in emission order rather than by variable name —
variable letters restart per chain and can repeat. A fully determined chain
contributes no equation at all.

`canonicalize` canonicalizes the molecular graph and updates every
atom reference in `Sg:`, `m:`, atom labels, and atom properties. Lipid trailer
tokens are retained and normalized with their referenced labels. Calling the
function twice returns the same string.

`smiles2name` always canonicalizes first. The reverse converter therefore does
not depend on the exact atom or branch order produced by `name2smiles`. Before
returning a name, it regenerates the structure and verifies canonical CXSMILES
equivalence.

## Depiction and toolkit behavior

CXSMILES support varies between chemistry toolkits. Use `smiles_for_depiction` before
rendering a string with `m:` blocks: it canonicalizes the complete string,
reindexes every CX field, and reduces each position variation to the two atoms
of the nearest unused side-chain single bond. Two endpoints are retained because
CDK ignores a one-endpoint multicenter bond and depicts its substituent as a
detached component.
The original string remains the analytical record. Renderers that ignore `Sg:`
see only the fixed scaffold and produce a truncated chain; use
`smiles_expand` or `name2structure` when a plain-SMILES
representative is required. That expansion does not turn an unknown position
into a measurement.

Round-tripping a stored string through a toolkit is lossy in a way worth
knowing about: the `|...|` block survives and is renumbered correctly, but the
trailer is a title and is dropped. A string that went in carrying
`constrain(a+b=15)` comes back with its `Sg:` markers intact and no length —
still well-formed, quietly less specific. Preserve the trailer yourself if a
toolkit rewrites the molecule.

## Limitations

- Carbohydrate sequences such as `Gal-Glc-Cer` do not determine glycosidic
  linkage positions and are not converted.
- A functional group on C1 of an acyl chain can change the linkage itself and
  is rejected.
- Confidence percentages never enter the structure. They are carried in the
  trailer as `dbPos`/`mPos` tokens and restored by `smiles2name`.
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

Regenerate the thematic conversion gallery with:

```bash
cargo test --test demo -- --ignored
```

This updates [`demo.html`](demo.html), including large live SVG depictions from
CDK Depict for each representable example.

## License

MIT

## References
- **Lipid shorthand and ID levels** — Liebisch G, et al. *Update on LIPID MAPS
  classification, nomenclature, and shorthand notation for MS-derived lipid
  structures.* J Lipid Res, 2020;61(12):1539–1555.
  [doi:10.1194/jlr.S120001025](https://doi.org/10.1194/jlr.S120001025) —
  Tables 1A/1B/1C are the modification and ring vocabulary of §1.3.
