# Lipid Shortname to/from (extended) CXSMILES conversion

Rationale: [Translating lipid shorthand notation into extended CXSMILES, and back](https://metabolomics.blog/2026/08/fair-lipid-representation-v2/)

`lipid_notation` converts [Shorthand2020](https://doi.org/10.1194/jlr.S120001025)
lipid shorthand names to SMILES/CXSMILES and back. It preserves structural
uncertainty instead of choosing positions or chain assignments that were not
measured.

It adheres to CXSMILES standards as closely as possible, and adds tokens after the `| |` pipes to carry information CXSMILES itself doesn't support or retain.

Check `demo.html` for examples.

## Quick start

```rust
use lipid_notation::{canonicalize, name2smiles, smiles2name, smiles_for_depiction};

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
as `;O2` are rejected unless the functional groups are specified. A chain
carrying both a *localized* group and *unlocalized* double bonds is rejected as
well, because the name does not say how the bonds distribute around the group; see
[Limitations](#limitations).

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
| geometry undetermined (`18:1(9)`) | yes: a bare `C=C` already means it | n/a |
| a double bond somewhere in a stretch (`18:1`) | no | partly: `Sg:` marks the stretch but has nowhere to record its length |
| a modification on one of many carbons (`;OH`) | no | yes: `m:` position-variation |
| chains known, sn assignment not (`16:0_18:1`) | no | badly: `RG:` over-generates and cannot coexist with `Sg:` |
| a weighted call (`Δ9 92%`) | no | no; it goes in the trailer instead |
| sum composition (`PC 34:1`) | no | no |

### The design rule

**Everything the standard covers goes inside the pipes, correct and complete.
Everything it does not goes after them, where no toolkit can mistake it for
structure.**

The part before and inside `|...|` is ordinary, conformant CXSMILES. A toolkit
that knows nothing about lipids parses it, renders it, and gets a chemically
valid molecule. Nothing lipid-specific is smuggled into a standard field, and no
standard field is given a private meaning.

### Standard CXSMILES blocks inside the pipes

| Encoding | Meaning |
|---|---|
| bare `C=C` | double-bond geometry undetermined: SMILES already means this, so nothing is added |
| `Sg:` | a double bond lies somewhere in this stretch |
| `m:` | this group attaches to one of these candidate atoms |
| `$snN$` atom labels | this atom is the sn-*N* attachment point |

### The trailer: everything CXSMILES cannot express

**This trailer format is a proposal, not a standard.** CXSMILES toolkits
tolerate an arbitrary string after the closing `|...|` pipes and typically
treat it as an opaque generic string, so this crate uses that space for a set
of `;`-separated tokens covering what CXSMILES cannot express: length
constraints on a `Sg:` run, regiochemistry left unknown because `Sg:` and
`RG:` cannot be combined, and confidence scores on bonds or atoms.

The trailer is a `;`-separated list of `name(argument)` tokens. It generalizes
the `a+b=15` size constraint proposed in [CDK CXSMILES](https://egonw.github.io/cdk-cxsmiles/templates.html#lipids-with-a-double-bond-somewhere-in-the-tail), into
something that can carry a second kind of statement:

| Token | Meaning |
|---|---|
| `constrain(a+b=15)` | the `Sg:` runs marked `a` and `b` span 15 carbons between them |
| `swappable(sn1,sn2)` | the chains at these labelled positions may be exchanged |
| `dbPos(sn1:9@92)` | the double bond at Δ9 on sn-1 was called with 92% confidence |
| `mPos(OH1:11OH@50,13OH@50)` | the group on the stub labelled `OH1` is at position 11 or 13, evenly split |

`dbPos` and `mPos` carry the bracketed consensus tail that instrument software
puts after a name, such as `FA 18:2(9,12) [DB sn1: Δ9 100%, Δ12 88%]`. No structure
format can hold a weighted call, so it is carried in the trailer instead —
metadata belongs where this crate's other metadata lives — and `smiles2name`
reconstructs the original tail from it.

Within a token, `,` separates positions and `|` separates *mutually exclusive*
candidates for one feature: `dbPos(sn1:5@100|14@50|15@50)` reads as "Δ5 for
certain, plus one more double bond that is either Δ14 or Δ15". Positions and
percentages are written `pos@percent` with no `Δ` and no spaces, because a
`.smi` reader splits the line on whitespace and would truncate the token; the
original spelling comes back on the way out.

These tokens only ever *refine* what the structure already says. An entry for a
position the SMILES commits to records how sure that call was. A set of `|`
alternatives corresponds to a double bond left inside an `Sg:` run, or a group
left on an `m:` stub. Narrowing "somewhere in this stretch" to "one of these,
with these odds" is more information, and still not a determination.

**These tokens are not official CXSMILES.** They sit in what a SMILES reader
treats as the *title* field, so a toolkit reads them as the molecule's name and
never as structure. CDK Depict prints them under the picture unless asked
not to, and RDKit stores them in `_Name` and drops them when it writes the
molecule out again. Strip the trailer and what remains is still valid CXSMILES
describing a chemically valid molecule; you lose the lipid semantics, not the
structure.

Two rules make the trailer safe to extend:

- **Tokens name things, never positions.** `constrain` names `Sg:` variables,
  and `swappable`, `dbPos` and `mPos` name `$...$` atom labels: a chain's
  `snN`, or an `m:` stub's own label. Toolkits maintain those labels across
  canonicalization: renumber the molecule and the label travels with its atom.
  An atom *index* in the trailer would silently rot the first time that
  happened, because nothing rewrites the title field.
- **Anything stated in the trailer is anchored in the pipes.** A `swappable`,
  `dbPos` or `mPos` token always comes with the `$...$` labels it names, so the
  `|...|` block is never empty while the trailer has something to say. That keeps the one-character
  check honest: **pipes mean something was undetermined.**

The same shape extends to anything else worth stating: the grammar is open, and
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
  (one label with *N* definitions substitutes each of its sites independently)
  `PC 16:0_18:1` also admitted `PC 16:0/16:0`, and a TG permitted 27 assignments
  where the name allows 6. Whether a given toolkit reads it that way or pairs
  definitions to sites positionally is itself the problem: nothing in the string
  says which is meant. Labelled positions plus `swappable` are unambiguous.
- **RDKit does not implement `RG:` at all**, not even the minimal case. Dropping
  it means every string this crate emits now parses in both toolkits.
- **`Sg:n` has nowhere to record a repeat count.** Hence `constrain(...)`.
- **`m:` candidate lists cannot see inside `Sg:` runs.** A chain with both an
  unlocalized modification and unlocalized double bonds can only offer the atoms
  physically present in the string, so the candidate list understates the
  ambiguity, which is the safer direction but still not the truth.
- **`ctu:` and `f:` look applicable and are not.** `ctu:` is a *query* feature
  for matching either geometry, and a bare `C=C` already says as much about a
  structure; `f:` groups components into one entity and expresses *and*, never
  the *or* that unresolved regiochemistry needs.
- **No construct expresses weighted alternatives**, so a 92% call and a 100%
  call are written identically in the structure. The percentages are not lost:
  they move to `dbPos`/`mPos` trailer tokens, which name the chain or the `m:`
  stub they refer to and never touch a CXSMILES field.


### Examples

Fully determined, so plain SMILES with no tail at all:

```
FA 18:1(9Z)  →  OC(=O)CCCCCCC/C=C\CCCCCCCC
```

One double bond, position not determined, so a flexible run plus its length
rather than a guessed position:

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

And the case that motivated the whole design: an unlocalized double bond *and*
an unresolved sn assignment, which the `RG:` encoding could not state together:

```
PC 16:0_18:1
→ C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CC=CC)COC(=O)CCCCCCCCCCCCCCC |$;;;;;;;;;;;;;sn2;;;;;;;;sn1$,Sg:n:16:a:ht,Sg:n:19:b:ht| swappable(sn1,sn2);constrain(a+b=15)
```

Worth noting: each chain gets its own equation, and equations are matched to
`Sg:` markers *positionally*, in emission order rather than by variable name.
Variable letters restart per chain and can repeat. A fully determined chain
contributes no equation at all.

`canonicalize` canonicalizes the molecular graph and updates every
atom reference in `Sg:`, `m:`, atom labels, and atom properties. Lipid trailer
tokens are retained and normalized with their referenced labels. Calling the
function twice returns the same string.

`smiles2name` always canonicalizes first. The reverse converter therefore does
not depend on the exact atom or branch order produced by `name2smiles`. Before
returning a name, it regenerates the structure and verifies canonical CXSMILES
equivalence.

#### Where name and string are not exactly interchangeable

A name that converts always yields a string that regenerates it, but the two
are not word-for-word equivalents. Round-tripping preserves the structure and
the ambiguity while normalizing the spelling to one form:

| Input | Recovered |
|---|---|
| `FA 18:0;ep(5)` | `FA 18:0;5Ep` |
| `FA 18:2;OH;OH` | `FA 18:2;(OH)2` |
| `CL 16:0/18:1(9)/18:2(9,12)/18:0` | `CL 18:2(9,12)/18:0/16:0/18:1(9)`; the two phosphatidyl halves are exchangeable across the central glycerol, so one ordering is chosen |
| `ST 27:1(5);3OH` | `ST`; the sterol form is chain-free, and the ring detail is implied by the class rather than carried in the name |

Compare canonical CXSMILES, not name strings, when testing equivalence.

Three further places where the string carries slightly more or less than the
name did:

**A `m:` candidate list understates ambiguity on a chain that also has `Sg:`
runs.** `FA 18:2;OH` emits `m:9:3.4.5.6.7.8`, six candidate carbons for a
hydroxyl that could sit on any of C2–C18. A `m:` list can only name atoms
physically present in the string, and the carbons hidden inside the flexible
runs are not. The error is in the safe direction (a reader concludes less than
the truth, never more), but it is still not the truth. The chain length in
`constrain(...)` is what tells you the list is partial.

**Confidence percentages never enter the structure.** `FA 18:2(9,12) [DB sn1:
Δ9 100%, Δ12 88%]` produces a fully committed SMILES plus
`dbPos(sn1:9@100,12@88)`. Strip the trailer, as any toolkit round-trip does,
and a 88% call is indistinguishable from a determination.

**`smiles2name` is not a general structure-to-name engine.** It recognizes
structures that match this crate's own templates. A valid lipid drawn some other
way, or a molecule outside the supported classes, returns `None` even when it is
chemically fine.

## Limitations

These names are refused: conversion returns `None` rather than inventing a
structure. The demo gallery (`cargo test --test demo -- --ignored`) renders the
same list. For names that *do* convert but whose string is not a word-for-word
equivalent of the name, see
[Where name and string are not exactly interchangeable](#where-name-and-string-are-not-exactly-interchangeable).

| Input | Why |
|---|---|
| `PC 34:1`, `TG 54:3` | Sum composition. Many chain combinations match, and picking one would be a fabrication. `class_needs_multi_chain` reports which classes this applies to. |
| `FA 18:1(9);O2` | A generic oxygen count names no functional group. Two oxygens could be two hydroxyls, a hydroperoxide, an epoxide plus a hydroxyl, or a ketone plus a hydroxyl. These are different molecules, not different positions of one molecule, so no `m:` block can cover them. Name the groups (`;(OH)2`, `;9Ep;OH`) and the oxygen count is no longer the obstacle. |
| `FA 18:0;1OMe`, `FA 18:1(9);1OH` | A group on C1 sits on the carboxyl carbon and redefines the linkage itself rather than decorating the chain. The result would no longer be the acyl species the class implies. |
| `Gal-Glc-Cer d18:1(4)/16:0` | A glycan sequence names the sugars and their order but not the glycosidic linkage positions or anomeric configurations. Those are bonds, not annotations, and there is nowhere in the string to leave them open. |
| `DG 16:0/18:1(9)/0:0` | `0:0` marks a vacant sn position in some dialects; this crate expects the vacancy to be expressed by the class (`DG` with two chains). |
| `FA 18:2;9OH`, `FA 18:2;9Ep`, `FA 19:1;[11-13cy3:0]` | A localized modification, epoxide, or ring combined with unlocalized double bonds. It's impossible to determine how many double bonds sit on which side of the fixed group. |

**A note on the last example: `FA 18:2;9OH`.** The molecule expressed by this shorthand name is a very unlikely result of an MS2 experiment. Assigning a `-OH` position implies that the number of C and H on both sides have been identified. Hence, the number of double bonds on both sides will be known. This could be easily expressed with `Sg:` blocks on both sides of the `-OH` fixed at C9. For the case of one C=C on both sides, the corresponding CXSMILES is
```
FA 18:2;9OH (one C=C is located on either sides of the -OH group)
→ OC(=O)CC=CCC(O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht,Sg:n:9:c:ht,Sg:n:12:d:ht| constrain(a+b=5);constrain(c+d=7)
```
This is a plausible results of a MS2 experiment. However, it can't be spelled out in shorthand notation. This is a (rare) limitation of the compact notation. Accordingly, `smiles2name` returns `None` for that string, since there is no name equivalent to it, and `name2smiles` returns `None` for `FA 18:2;9OH`, rather than emitting one arbitrary distribution of the two double bonds.

The rejection triggers only when the ambiguity is real: when there is room for
a double bond in front of the fixed group as well as behind it. A group close
enough to C1 leaves no such room, so exactly one distribution is consistent with
the name and the string is exact: `FA 18:2;3OH` converts. So do names that
resolve the question by localizing the double bonds (`FA 18:2(11,14);9OH`), the
same group on a saturated chain (`FA 18:0;9OH`), and unlocalized double bonds on
an unmodified chain (`FA 18:2`).

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
`constrain(a+b=15)` comes back with its `Sg:` markers intact and no length:
still well-formed, quietly less specific. Preserve the trailer yourself if a
toolkit rewrites the molecule.

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
`cargo test --test testdata -- --ignored` and review the corpus diff. For external
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
- **Lipid shorthand and ID levels**: Liebisch G, et al. *Update on LIPID MAPS
  classification, nomenclature, and shorthand notation for MS-derived lipid
  structures.* J Lipid Res, 2020;61(12):1539–1555.
  [doi:10.1194/jlr.S120001025](https://doi.org/10.1194/jlr.S120001025).
  Tables 1A/1B/1C are the modification and ring vocabulary of §1.3.
- **Blog post**: [Translating lipid shorthand notation into extended CXSMILES, and back](https://metabolomics.blog/2026/08/fair-lipid-representation-v2/)
