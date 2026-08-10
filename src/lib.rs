//! Convert Shorthand2020 lipid names to SMILES/CXSMILES and back.
//!
//! Plain SMILES represents fully determined structures. CXSMILES fields retain
//! lipid ambiguity without inventing positions:
//!
//! - `Sg:` marks an unlocalized double-bond region;
//! - `m:` lists candidate atoms for an unlocalized modification;
//! - `$snN$` labels plus `swappable(...)` retain unresolved sn assignment;
//! - `constrain(...)` records the length of an `Sg:` region.
//!
//! ```
//! use lipid_notation::name2smiles;
//!
//! assert_eq!(
//!     name2smiles("FA 18:1(9Z)").as_deref(),
//!     Some(r"OC(=O)CCCCCCC/C=C\CCCCCCCC")
//! );
//!
//! // The double-bond position is unknown, so the result carries an Sg region.
//! assert_eq!(
//!     name2smiles("FA 18:1").as_deref(),
//!     Some("OC(=O)CC=CC |Sg:n:3:a:ht,Sg:n:6:b:ht| constrain(a+b=15)")
//! );
//! ```
//!
//! [`smiles2name`] accepts equivalent atom and branch orders by canonicalizing
//! before recognition. [`expand_cxsmiles_for_depiction`] turns variable `Sg:`
//! regions into one concrete plain-SMILES representative.
//!
//! [CXSMILES]: https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html

mod cxsmiles;
mod forward;
mod nomenclature;
mod reverse;

pub use cxsmiles::canonicalize_cxsmiles;
pub use forward::{expand_cxsmiles_for_depiction, ChainAtoms, LipidStructure};

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
    forward::generate_smiles(name)
}

/// Reads a SMILES/CXSMILES string back into a Shorthand2020 lipid name.
///
/// This inverts [`name2smiles`]; it is not a general SMILES parser, and
/// returns `None` for a structure this crate would not have written.
///
/// The input is canonicalized first, so equivalent SMILES atom/branch orders
/// are accepted. **Every answer is proved before it is returned:** the name is
/// fed back through [`name2smiles`] and must regenerate the same canonical
/// CXSMILES, so this can lose coverage but cannot return a different structure.
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
/// * accepted aliases normalize — `;ep(5)` returns as `;5Ep`.
pub fn smiles2name(smiles: &str) -> Option<String> {
    reverse::parse_smiles(smiles)
}

/// [`name2smiles`] plus, for every chain, the atom index of each of its
/// carbons — the mapping a UI needs to highlight which part of a structure a
/// given MS2 fragment came from.
///
/// Differs from [`name2smiles`] in two depiction-driven ways: the `Sg:`
/// markers are already expanded to one representative even split, and
/// `_`-joined names are built as one representative assignment and flagged via
/// [`LipidStructure::regio_resolved`].
pub fn name2structure(name: &str) -> Option<LipidStructure> {
    forward::generate_structure(name)
}

/// Whether a lipid class inherently spans several acyl/alkyl chains, so a
/// bare sum-composition token can't be resolved to one real structure.
///
/// Exposed because callers that localize double bonds need the same test:
/// localizing a position within a *sum* composition would treat two real
/// chains as if they were one.
pub fn class_needs_multi_chain(class: &str) -> bool {
    forward::class_needs_multi_chain(class)
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
