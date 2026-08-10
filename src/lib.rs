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
//! use shorthand2smiles::name2smiles;
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
//!     Some("OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| a+b=15")
//! );
//! ```
//!
//! The presence of a `|...|` tail is itself the signal: **pipes mean
//! something was undetermined.**
//!
//! # The two directions
//!
//! * [`name2smiles`] — implemented, and the reason this crate exists.
//! * [`smiles2name`] — **not implemented yet**; always returns `None`. See its
//!   docs for why it is worth building.
//!
//! # Which block means what
//!
//! | block | says | emitted when |
//! |---|---|---|
//! | *(none)* | geometry undetermined | always — a bare `C=C` already means "cis or trans, not determined" |
//! | `Sg:` | "the double bond is somewhere in this stretch" | a chain declares more double bonds than it localizes |
//! | `m:` | "this group attaches somewhere on this chain" | a modification is written with no position (`;OH`) |
//! | `RG:` | "either chain could be at either position" | chains joined with `_`, and none of them needs `Sg:`/`m:` |
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
//! ignores `Sg:` (handing back a molecule missing most of its chain) and
//! rejects `RG:` outright. Use [`expand_cxsmiles_for_depiction`] to get a
//! plain SMILES first, or [`name2structure`] if you also need to know which
//! atoms belong to which chain.
//!
//! [CXSMILES]: https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html

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
/// use shorthand2smiles::name2smiles;
///
/// // sn-position unknown -> the chains become R-group alternatives
/// let s = name2smiles("DG 16:0_18:1(9)").unwrap();
/// assert!(s.contains("RG:_R1="));
///
/// // unsupported: sum composition only
/// assert_eq!(name2smiles("PC 34:1"), None);
/// ```
pub fn name2smiles(name: &str) -> Option<String> {
    smiles::lipid_name_to_smiles(name)
}

/// Parses a SMILES/CXSMILES string back into a Shorthand2020 lipid name.
///
/// **Not implemented — always returns `None`.** The signature is fixed and
/// the tests are written (see `smiles2name_round_trip`, currently
/// `#[ignore]`d), so filling this in should not disturb any caller.
///
/// It is worth building for two independent reasons:
///
/// 1. **Validation.** The forward mapping has no parser-level check that a
///    string means what the name meant. Round-tripping every name through
///    SMILES and back would catch index drift, miscounted chains and
///    malformed blocks in one pass — the class of bug that literal-comparison
///    tests structurally cannot catch, and that
///    [`name2smiles`] shipped three of before review caught them.
/// 2. **Ingestion.** Reading structures produced by other tools is currently
///    impossible, which rules out using this crate to normalize a
///    third-party lipid library.
///
/// The hard part is not SMILES parsing but *recovering the ambiguity*:
/// mapping `Sg:` runs back to "N declared, K given", and an `m:` candidate
/// list back to a position-less modification token, without silently
/// promoting a guess to a determination.
pub fn smiles2name(_smiles: &str) -> Option<String> {
    None
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

    #[test]
    fn smiles2name_is_not_implemented_yet() {
        // Documents the current state rather than the intended one. Delete
        // this test — don't edit it — when `smiles2name` lands, and remove
        // the `#[ignore]` from `smiles2name_round_trip` below.
        assert_eq!(smiles2name("OC(=O)CCCCCCCCCCCCCCCCC"), None);
    }

    /// The check `smiles2name` exists to make possible: every name that
    /// produces a structure must be recoverable from that structure.
    #[test]
    #[ignore = "smiles2name not implemented yet"]
    fn smiles2name_round_trip() {
        for name in [
            "FA 18:0",
            "FA 18:1(9Z)",
            "FA 18:1(9)",
            "FA 18:1",
            "FA 18:0;OH",
            "FA 18:0;5OH",
            "PC 16:0/18:1(9)",
            "PC 16:0_18:1(9)",
            "TG 16:0/18:1(9)/18:2(9,12)",
            "Cer d18:1(4)/16:0",
        ] {
            let smi = name2smiles(name).unwrap_or_else(|| panic!("{name} should resolve"));
            assert_eq!(
                smiles2name(&smi).as_deref(),
                Some(name),
                "{name} did not survive the round trip via {smi}"
            );
        }
    }
}
