//! Shorthand2020 lipid names ↔ SMILES/CXSMILES.
//!
//! Lipidomics reports structure at whatever level the evidence supports:
//! `PC 34:1` concedes everything but the sum composition, `PC 16:0_18:1`
//! names the chains but not their positions, `PC 16:0/18:1(9Z)` concedes
//! nothing. SMILES only speaks the last of those. This crate bridges the gap
//! by emitting [CXSMILES] blocks that state precisely what was *not*
//! determined, rather than guessing a fully specified structure and
//! presenting it as fact.
//!
//! ```
//! use lipid_notation::name2smiles;
//!
//! // Fully determined -> plain, ordinary SMILES.
//! assert_eq!(
//!     name2smiles("FA 18:1(9Z)").as_deref(),
//!     Some(r"OC(=O)CCCCCCC/C=C\CCCCCCCC")
//! );
//!
//! // Double bond present, position unknown -> a variable-length run plus a
//! // size constraint, never a guessed position.
//! assert_eq!(
//!     name2smiles("FA 18:1").as_deref(),
//!     Some("OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)")
//! );
//! ```
//!
//! The presence of a `|...|` tail is itself the signal: **pipes mean
//! something was undetermined.** The text after the closing pipe is a
//! `;`-separated list of tokens that are this crate's own, not CXSMILES —
//! `constrain(a+b=15)` sizes a flexible run, `swappable(sn1,sn2)` says an sn
//! assignment is one arbitrary choice. See `dev/extension.md`.
//!
//! [`smiles2name`] reads that back. The hard part is not SMILES parsing but
//! *recovering the ambiguity* — mapping `Sg:` runs back to "N declared, K
//! given" without silently promoting a guess to a determination — so every
//! answer it gives is re-generated and compared before being returned.
//!
//! # Which block means what
//!
//! | block | says | emitted when |
//! |---|---|---|
//! | *(none)* | geometry undetermined | always — a bare `C=C` already means "cis or trans, not determined" |
//! | `Sg:` | "the double bond is somewhere in this stretch" | a chain declares more double bonds than it localizes |
//! | `m:` | "this group attaches somewhere on this chain" | a modification is written with no position (`;OH`) |
//! | `$snN$` + `swappable(...)` | "either chain could be at either position" | chains joined with `_` |
//!
//! Two CXSMILES constructs are deliberately **not** used. `ctu:` is a query
//! feature for matching either configuration when searching — a plain `C=C`
//! already says the same thing about a structure. `f:` groups components into
//! one entity (a salt, a hydrate) and expresses *and*, never *or*, so it
//! cannot say "this chain **or** that chain sits at sn-1". Earlier revisions
//! used both; `dev/SMILES.md` §3 records what went wrong.
//!
//! # Depiction
//!
//! Stored strings are rigorous, not drawable by every tool: RDKit silently
//! ignores `Sg:`, handing back a molecule missing most of its chain. Use
//! [`expand_cxsmiles_for_depiction`] to get a
//! plain SMILES first, or [`name2structure`] if you also need to know which
//! atoms belong to which chain.
//!
//! [CXSMILES]: https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html

mod from_smiles;
mod nomenclature;
mod smiles;

pub use smiles::{expand_cxsmiles_for_depiction, ChainAtoms, LipidStructure};

/// Converts a Shorthand2020 lipid name into SMILES, or CXSMILES when
/// something about the structure was not determined.
///
/// Returns `None` when the name cannot be turned into a structure honestly:
///
/// * the headgroup has no template here (unsupported class);
/// * multi-chain **shorthand** (`PC 34:1`, `TG 54:3`) — the sum composition
///   has many chain realizations and picking one would be a fabrication.
///   Supply explicit chains instead (`PC 16:0/18:1`, `PC 16:0_18:1`);
/// * a chain's oxygen count is declared generically (`;O2`) with no
///   `;OH`/`;oxo`/`;COOH` breakdown, so there is no position hypothesis at
///   all.
///
/// Any bracketed confidence tail (`FA 18:2 [DB sn1: Δ9 92%]`) is stripped
/// first; no structure format can carry a weighted distribution.
///
/// ```
/// use lipid_notation::name2smiles;
///
/// // sn-position unknown -> the linking atoms are labelled and the trailing
/// // token says their assignment may be permuted
/// let s = name2smiles("DG 16:0_18:1(9)").unwrap();
/// assert!(s.contains("swappable(sn1,sn2)"));
///
/// // unsupported: sum composition only
/// assert_eq!(name2smiles("PC 34:1"), None);
/// ```
pub fn name2smiles(name: &str) -> Option<String> {
    smiles::lipid_name_to_smiles(name)
}

/// Reads a SMILES/CXSMILES string back into a Shorthand2020 lipid name.
///
/// This inverts [`name2smiles`]; it is not a general SMILES parser, and
/// returns `None` for a structure this crate would not have written.
///
/// **Every answer is proved before it is returned.** The name is fed back
/// through [`name2smiles`] and must regenerate the input string exactly, so
/// this can lose coverage but cannot hand back a name meaning something
/// other than the structure given.
///
/// ```
/// use lipid_notation::{name2smiles, smiles2name};
///
/// let smi = name2smiles("FA 18:1(9Z)").unwrap();
/// assert_eq!(smiles2name(&smi).as_deref(), Some("FA 18:1(9Z)"));
///
/// // Undetermined positions survive the trip: the Sg: run comes back as a
/// // double bond count with no position, exactly as it went in.
/// let smi = name2smiles("FA 18:1").unwrap();
/// assert_eq!(smiles2name(&smi).as_deref(), Some("FA 18:1"));
///
/// // Not something this crate wrote.
/// assert_eq!(smiles2name("c1ccccc1"), None);
/// ```
///
/// # Where the round trip is lossy
///
/// [`name2smiles`] is not injective, so the name returned is not always the
/// name you started from — but it always regenerates the same string:
///
/// * a trailing empty slot is dropped, so `DG 16:0/0:0` and `MG 16:0` are the
///   same string and both return `MG 16:0`;
/// * older spellings canonicalize — `;ep(5)` returns as `;5Ep`.
pub fn smiles2name(smiles: &str) -> Option<String> {
    from_smiles::smiles_to_name(smiles)
}

/// [`name2smiles`] plus, for every chain, the atom index of each of its
/// carbons — the mapping a UI needs to highlight which part of a structure a
/// given MS2 fragment came from.
///
/// Differs from [`name2smiles`] in two depiction-driven ways: the `Sg:`
/// markers are already expanded to one representative even split, and
/// `_`-joined names are built as if they were `/`-joined (flagged via
/// [`LipidStructure::regio_resolved`], since a Markush scheme has no concrete
/// atoms to point at). Both trades are documented in `dev/SMILES.md` §2.3.
pub fn name2structure(name: &str) -> Option<LipidStructure> {
    smiles::lipid_name_to_structure(name)
}

/// Whether a lipid class inherently spans several acyl/alkyl chains, so a
/// bare sum-composition token can't be resolved to one real structure.
///
/// Exposed because callers that localize double bonds need the same test:
/// localizing a position within a *sum* composition would treat two real
/// chains as if they were one.
pub fn class_needs_multi_chain(class: &str) -> bool {
    smiles::class_needs_multi_chain(class)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name2smiles_is_the_documented_entry_point() {
        assert_eq!(
            name2smiles("FA 18:0").as_deref(),
            Some("OC(=O)CCCCCCCCCCCCCCCCC")
        );
    }
}
