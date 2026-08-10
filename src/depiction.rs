//! CXSMILES adjustments for predictable lipid depictions.

use std::collections::{HashSet, VecDeque};

use chematic_core::{AtomIdx, BondOrder, Element, Molecule};
use chematic_smiles::parse;

use crate::cxsmiles::bracket_bare_wildcards;
use crate::reverse::{reachable_from, single_bond};

/// Prepares a lipid CXSMILES string for a deterministic depiction.
///
/// An `m:` position-variation block correctly says that a disconnected
/// substituent may attach at any of several atoms. Depiction tools must still
/// choose one attachment to draw, and their choice can be arbitrary (and may
/// put a group on a ring when a side-chain position is more helpful). This
/// function first canonicalizes the complete CXSMILES, reindexing every CX
/// field, and then replaces each `m:` candidate list with the two endpoints of
/// one representative single bond. The bond nearest the chain's carbonyl or
/// other headgroup attachment is preferred, and multiple unlocalized
/// modifications use different bonds.
///
/// `Sg:`, atom labels, and lipid trailer tokens are retained. The returned
/// CXSMILES is therefore a drawing representative, not a replacement for the
/// original analytical record. Malformed or non-CXSMILES input is returned
/// unchanged.
pub fn smiles_for_depiction(smi: &str) -> String {
    let Some(canonical) = crate::canonicalize(smi) else {
        return smi.to_string();
    };
    let Some((base, rest)) = canonical.split_once(" |") else {
        return canonical;
    };
    let Some((fields, trailer)) = rest.split_once('|') else {
        return canonical;
    };

    let molecule = parse(&bracket_bare_wildcards(base)).ok();
    let mut selected = HashSet::new();
    let fields = fields
        .split(',')
        .map(|field| rewrite_position_variation(field, molecule.as_ref(), &mut selected))
        .collect::<Vec<_>>()
        .join(",");
    format!("{base} |{fields}|{trailer}")
}

fn rewrite_position_variation(
    field: &str,
    molecule: Option<&Molecule>,
    selected: &mut HashSet<usize>,
) -> String {
    let Some(rest) = field.strip_prefix("m:") else {
        return field.to_string();
    };
    let Some((floating, candidates)) = rest.split_once(':') else {
        return field.to_string();
    };
    let Some(floating) = floating.parse::<usize>().ok() else {
        return field.to_string();
    };
    let candidates = candidates
        .split('.')
        .map(str::parse::<usize>)
        .collect::<Result<Vec<_>, _>>();
    let Ok(candidates) = candidates else {
        return field.to_string();
    };
    let Some((first, second)) = preferred_bond(molecule, floating, &candidates, selected) else {
        return field.to_string();
    };

    selected.extend([first, second]);
    format!("m:{floating}:{first}.{second}")
}

fn preferred_bond(
    molecule: Option<&Molecule>,
    floating: usize,
    candidates: &[usize],
    selected: &HashSet<usize>,
) -> Option<(usize, usize)> {
    let molecule = molecule?;
    if floating >= molecule.atom_count()
        || candidates
            .iter()
            .any(|candidate| *candidate >= molecule.atom_count())
    {
        return None;
    }

    let available = candidates
        .iter()
        .copied()
        .filter(|candidate| !selected.contains(candidate))
        .collect::<HashSet<_>>();

    let excluded = reachable_from(molecule, floating, &HashSet::new());
    let carbonyls = (0..molecule.atom_count())
        .filter(|atom| !excluded.contains(atom) && is_carbonyl_carbon(molecule, *atom))
        .collect::<HashSet<_>>();
    let headgroup_atoms = (0..molecule.atom_count())
        .filter(|atom| !excluded.contains(atom) && is_headgroup_atom(molecule, *atom))
        .collect::<HashSet<_>>();

    let mut bonds = Vec::new();
    for &first in &available {
        if excluded.contains(&first) {
            continue;
        }
        for (neighbor, bond) in molecule.neighbors(AtomIdx(first as u32)) {
            let second = neighbor.0 as usize;
            if first >= second
                || !available.contains(&second)
                || !single_bond(molecule.bond(bond).order)
            {
                continue;
            }
            let first_score = site_score(
                molecule,
                first,
                candidates,
                &carbonyls,
                &headgroup_atoms,
                &excluded,
            );
            let second_score = site_score(
                molecule,
                second,
                candidates,
                &carbonyls,
                &headgroup_atoms,
                &excluded,
            );
            let (near, far, score) = if first_score <= second_score {
                (first, second, first_score)
            } else {
                (second, first, second_score)
            };
            bonds.push((score, near, far));
        }
    }
    bonds.sort_unstable();
    bonds.first().map(|(_, near, far)| (*near, *far))
}

fn site_score(
    molecule: &Molecule,
    site: usize,
    candidates: &[usize],
    carbonyls: &HashSet<usize>,
    headgroup_atoms: &HashSet<usize>,
    excluded: &HashSet<usize>,
) -> (usize, usize, usize) {
    (
        graph_distance(molecule, site, carbonyls, excluded).unwrap_or(usize::MAX),
        graph_distance(molecule, site, headgroup_atoms, excluded).unwrap_or(usize::MAX),
        candidates
            .iter()
            .position(|candidate| *candidate == site)
            .unwrap_or(usize::MAX),
    )
}

fn is_carbonyl_carbon(molecule: &Molecule, atom: usize) -> bool {
    molecule.atom(AtomIdx(atom as u32)).element == Element::C
        && molecule
            .neighbors(AtomIdx(atom as u32))
            .any(|(neighbor, bond)| {
                molecule.atom(neighbor).element == Element::O
                    && molecule.bond(bond).order == BondOrder::Double
            })
}

fn is_headgroup_atom(molecule: &Molecule, atom: usize) -> bool {
    let element = molecule.atom(AtomIdx(atom as u32)).element;
    element != Element::C && molecule.neighbors(AtomIdx(atom as u32)).count() > 1
}

fn graph_distance(
    molecule: &Molecule,
    start: usize,
    targets: &HashSet<usize>,
    excluded: &HashSet<usize>,
) -> Option<usize> {
    if targets.is_empty() || excluded.contains(&start) {
        return None;
    }
    let mut seen = HashSet::from([start]);
    let mut todo = VecDeque::from([(start, 0)]);
    while let Some((atom, distance)) = todo.pop_front() {
        if targets.contains(&atom) {
            return Some(distance);
        }
        for (neighbor, _) in molecule.neighbors(AtomIdx(atom as u32)) {
            let neighbor = neighbor.0 as usize;
            if !excluded.contains(&neighbor) && seen.insert(neighbor) {
                todo.push_back((neighbor, distance + 1));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depicts_an_unlocalized_hydroxyl_next_to_the_fatty_acid_headgroup() {
        let cxsmiles = crate::name2smiles("FA 20:2(5,8);[11-15cy5;13OH];OH").unwrap();
        let depiction = smiles_for_depiction(&cxsmiles);
        assert!(depiction.contains("m:1:17.16|"), "{depiction}");
    }

    #[test]
    fn finds_the_headgroup_site_after_canonicalization() {
        let cxsmiles = crate::name2smiles("FA 20:2(5,8);[11-15cy5;13OH];OH").unwrap();
        let canonical = crate::canonicalize(&cxsmiles).unwrap();
        let depiction = smiles_for_depiction(&canonical);
        assert!(depiction.contains("m:1:17.16|"), "{depiction}");
        assert_eq!(
            depiction,
            smiles_for_depiction(&cxsmiles),
            "native and canonical inputs must produce the same depiction form"
        );
    }

    #[test]
    fn keeps_non_position_variation_cx_fields() {
        let cxsmiles = crate::name2smiles("PC 16:0_18:1").unwrap();
        assert_eq!(
            smiles_for_depiction(&cxsmiles),
            crate::canonicalize(&cxsmiles).unwrap()
        );
    }
}
