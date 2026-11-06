# Architecture

`lipidoracle-smiles` (`lipid_notation`) is a single-crate Rust library translating between Shorthand2020 lipid names and ambiguity-preserving SMILES/CXSMILES strings powered by `chematic-core` and `chematic-smiles`.

## Module Boundaries

- **`src/lib.rs`**: Public crate root exposing `name2smiles`, `smiles2name`, `name2structure`, `class_needs_multi_chain`, `canonicalize`, `smiles_for_depiction`, `smiles_expand`, `LipidStructure`, and `ChainAtoms`.
- **`src/nomenclature.rs`**: Lexer and parser for Shorthand2020 lipid class prefixes, headgroup stereochemistry, sn-separators (`/` vs. `_`), acyl/alkyl/alkenyl chain descriptors, double bonds, and functional-group modifications.
- **`src/forward.rs`**: Synthesizes molecular graphs and CXSMILES annotations (`Sg:`, `m:`, `$snN$`, `constrain(...)`, `swappable(...)`) from parsed lipid descriptors, and builds atom-index maps (`LipidStructure`, `ChainAtoms`) for MS2 fragment highlighting.
- **`src/reverse.rs`**: Recognizes lipid headgroup templates, chain attachments, double bonds, and CXSMILES ambiguity markers from an arbitrary SMILES/CXSMILES input and proves correctness by verifying that `name2smiles` reproduces the exact canonical CXSMILES.
- **`src/cxsmiles.rs`**: Parses, canonicalizes, and re-indexes CXSMILES atom labels, polymer repeat units (`Sg:`), positional variation bonds (`m:`), and constraint/swappable trailers across `chematic-smiles` canonicalization.
- **`src/depiction.rs`**: Converts `m:` positional-variation CXSMILES into deterministic 2D-depictable structures (`smiles_for_depiction`).
- **`src/consensus.rs`**: Serializes and deserializes bracketed localization-probability tails (`[DB ...]`) into `dbPos`/`mPos` CXSMILES trailer tokens so consensus calls survive round-tripping.

## Pipeline Stages

### Forward (`name2smiles` / `name2structure`)

1. **Consensus Extraction**: Strip and parse any trailing bracketed consensus annotation (`src/consensus.rs`).
2. **Shorthand2020 Parsing**: Parse lipid class, headgroup variant, chain separator (`/` resolved vs. `_` unresolved), and per-chain carbon/double-bond/modification specs (`src/nomenclature.rs`). Reject sum-only multi-chain compositions (`PC 34:1`) that cannot be represented without inventing chains.
3. **Template & Chain Assembly**: Instantiate the headgroup SMILES template, graft each chain with exact stereo (`@`/`@@`, `/` `\`) or `Sg:`/`m:` uncertainty regions, and attach `constrain(...)`, `swappable(...)`, and `dbPos`/`mPos` trailer tokens (`src/forward.rs`).

### Reverse (`smiles2name`)

1. **CXSMILES Canonicalization**: Parse the SMILES graph and CXSMILES trailer and canonicalize atom numbering via `cxsmiles::canonicalize`.
2. **Subgraph Recognition**: Match headgroup core, trace attached acyl/alkyl chains, and decode double-bond geometry, ring/oxygen bridges, and `Sg:`/`m:`/`swappable` ambiguity layers into a candidate Shorthand2020 name (`src/reverse.rs`).
3. **Round-Trip Verification**: Run the candidate name through `name2smiles` and `cxsmiles::canonicalize`; return `Some(name)` only if the regenerated canonical CXSMILES matches the input identically.
