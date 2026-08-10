//! Reverse conversion from SMILES/CXSMILES to Shorthand2020 lipid notation.
//!
//! Headgroups and chains are recognized from the molecule graph, so any atom
//! or branch order works — this crate's own output, a canonicalized form, or a
//! third-party spelling of the same structure. Every candidate name is
//! regenerated and compared in canonical CXSMILES form before it is returned,
//! so a bug here costs coverage but cannot return a name that means something
//! other than the structure given.

use crate::canonicalize;
use crate::forward::{generate_smiles, trailer_equations, trailer_is_swappable, SUBSTITUENTS};
use chematic_core::{AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_smiles::{canonical_smiles, parse};
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

/// A structure read back far enough to name one chain.
#[derive(Debug, Default, Clone)]
struct Chain {
    /// `""`, `"O-"`, `"P-"`, `"d"` or `"t"`.
    prefix: &'static str,
    carbon: u32,
    /// Localized double bonds: Δ-position and geometry, if the string
    /// carried one.
    db: Vec<(u32, Option<char>)>,
    /// Double bonds the `Sg:` scaffold stands in for — counted, never placed.
    unlocalized_db: usize,
    /// `(position, abbreviation)`, position `0` meaning an `m:` block said
    /// "present, not determined".
    mods: Vec<(u32, &'static str)>,
    /// `(first carbon, last carbon, is_epoxide)`.
    rings: Vec<(u32, u32, bool)>,
}

/// The CXSMILES tail, split into the parts that mean something here: the
/// standard `|...|` blocks and the lipid-specific trailing token list.
struct Cx<'a> {
    base: &'a str,
    /// Atom index of each `Sg:n:` marker, in emission order.
    sg: Vec<usize>,
    /// `(floating atom index, candidate atom indices)` for each `m:` block.
    m: Vec<(usize, Vec<usize>)>,
    /// Whether a `swappable(...)` token says the sn assignment shown is one
    /// arbitrary choice — i.e. whether the name joined its chains with `_`.
    swappable: bool,
    /// `(term count, sum)` per `constrain(...)` token, in order. They pair
    /// with `sg` positionally — see `smiles_expand`.
    equations: Vec<(usize, usize)>,
}

impl<'a> Cx<'a> {
    fn split(smi: &'a str) -> Cx<'a> {
        let Some((base, tail)) = smi.split_once(" |") else {
            return Cx {
                base: smi,
                sg: Vec::new(),
                m: Vec::new(),
                swappable: false,
                equations: Vec::new(),
            };
        };
        let (blocks, trailer) = match tail.split_once('|') {
            Some((b, c)) => (b, c),
            None => (tail, ""),
        };

        let mut sg = Vec::new();
        let mut m = Vec::new();
        for block in blocks.split(',') {
            if let Some(rest) = block.strip_prefix("Sg:n:") {
                if let Some(Ok(atom)) = rest.split(':').next().map(str::parse) {
                    sg.push(atom);
                }
            } else if let Some(rest) = block.strip_prefix("m:") {
                if let Some((idx, sites)) = rest.split_once(':') {
                    if let Ok(idx) = idx.parse() {
                        m.push((
                            idx,
                            sites.split('.').filter_map(|s| s.parse().ok()).collect(),
                        ));
                    }
                }
            }
        }

        // The trailer grammar lives with the emitter, in `smiles.rs`.
        let swappable = trailer_is_swappable(trailer);
        let equations = trailer_equations(trailer);

        Cx {
            base,
            sg,
            m,
            swappable,
            equations,
        }
    }
}

/// Writes a chain back out as its shorthand token, e.g. `18:1(9Z);5OH`.
fn format_chain(chain: &Chain) -> String {
    // A ring's double bonds are declared inside the ring and counted there,
    // so they belong to neither the chain's `C:DBE` count nor its position
    // list. A plasmenyl chain's C1=C2 is likewise implied by the `P-` prefix.
    let in_ring = |pos: u32| {
        chain
            .rings
            .iter()
            .any(|(s, e, bridged)| !*bridged && pos >= *s && pos < *e)
    };
    let placed: Vec<&(u32, Option<char>)> = chain
        .db
        .iter()
        .filter(|(p, _)| !(in_ring(*p) || (chain.prefix == "P-" && *p == 1)))
        .collect();

    let mut out = format!(
        "{}{}:{}",
        chain.prefix,
        chain.carbon,
        placed.len() + chain.unlocalized_db
    );
    if !placed.is_empty() {
        let list: Vec<String> = placed
            .iter()
            .map(|(p, g)| match g {
                Some(g) => format!("{p}{g}"),
                None => p.to_string(),
            })
            .collect();
        out.push_str(&format!("({})", list.join(",")));
    }

    // Rings come before the other groups, with any group inside their span
    // written within their brackets.
    for &(start, end, bridged) in &chain.rings {
        if bridged {
            out.push_str(&format!(";{start}Ep"));
            continue;
        }
        let ring_db: Vec<String> = chain
            .db
            .iter()
            .filter(|(p, _)| *p >= start && *p < end)
            .map(|(p, g)| match g {
                Some(g) => format!("{p}{g}"),
                None => p.to_string(),
            })
            .collect();
        let mut inner = format!("{start}-{end}cy{}:{}", end - start + 1, ring_db.len());
        if !ring_db.is_empty() {
            inner.push_str(&format!("({})", ring_db.join(",")));
        }
        for (pos, abbr) in chain.mods.iter().filter(|(p, _)| *p >= start && *p <= end) {
            inner.push_str(&format!(";{pos}{abbr}"));
        }
        out.push_str(&format!(";[{inner}]"));
    }

    // Then the rest, grouped by abbreviation in the table's own order so a
    // name comes back out the way the paper writes it.
    for (abbr, _) in SUBSTITUENTS {
        let positioned: Vec<String> = chain
            .mods
            .iter()
            .filter(|(p, a)| a == abbr && *p != 0 && !in_ring_span(chain, *p))
            .map(|(p, a)| format!("{p}{a}"))
            .collect();
        if !positioned.is_empty() {
            out.push_str(&format!(";{}", positioned.join(",")));
        }
        let loose = chain
            .mods
            .iter()
            .filter(|(p, a)| a == abbr && *p == 0)
            .count();
        if loose > 0 {
            // The paper parenthesizes a group when it occurs more than once,
            // and whenever its abbreviation carries a digit of its own.
            if loose > 1 || abbr.contains(|c: char| c.is_ascii_digit()) {
                out.push_str(&format!(";({abbr})"));
                if loose > 1 {
                    out.push_str(&loose.to_string());
                }
            } else {
                out.push_str(&format!(";{abbr}"));
            }
        }
    }
    out
}

fn in_ring_span(chain: &Chain, pos: u32) -> bool {
    chain
        .rings
        .iter()
        .any(|(s, e, bridged)| !*bridged && pos >= *s && pos <= *e)
}

/// `(class, headgroup tail, carries an N-acyl chain)`.
const SPHINGO: &[(&str, &str, bool)] = &[
    ("SM", "OP(=O)([O-])OCC[N+](C)(C)C", true),
    ("HexCer", "OC1OC(CO)C(O)C(O)C1O", true),
    ("IPC", "OP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", true),
    ("CerP", "OP(=O)(O)O", true),
    ("Cer", "O", true),
    ("S1P", "OP(=O)(O)O", false),
    ("Sph", "O", false),
];

/// Assembles the chain tokens into sn order, writing `0:0` for a gap that
/// some later position fills. The separator is `_` when a `swappable` token
/// said the sn assignment shown is arbitrary, `/` when it is determined.
fn join_slots(mut slots: Vec<(usize, Option<String>)>, sep: &str) -> Option<(usize, String)> {
    slots.sort_by_key(|(sn, _)| *sn);
    let last = slots.iter().rposition(|(_, c)| c.is_some())?;
    let tokens: Vec<String> = slots[..=last]
        .iter()
        .map(|(_, c)| c.clone().unwrap_or_else(|| "0:0".to_string()))
        .collect();
    Some((last + 1, tokens.join(sep)))
}

/// Reads a SMILES/CXSMILES string back into a Shorthand2020 name, or `None`
/// when no name this generator would accept regenerates it.
pub(crate) fn parse_smiles(smi: &str) -> Option<String> {
    let smi = smi.trim();
    let canonical = canonicalize(smi)?;

    let matches = |candidate: &String| {
        generate_smiles(candidate)
            .and_then(|generated| canonicalize(&generated))
            .as_deref()
            == Some(canonical.as_str())
    };

    // Recognition works on the molecule graph, so the input's atom and
    // branch order does not matter. Every hypothesis is verified in
    // canonical space before it is returned.
    graph_candidates(&canonical).into_iter().find(matches)
}
struct GraphTemplate {
    class: &'static str,
    chains: usize,
    sn: Vec<usize>,
    head_signature: String,
}

fn graph_candidates(canonical: &str) -> Vec<String> {
    let cx = Cx::split(canonical);
    let Ok(mol) = parse(cx.base) else {
        return Vec::new();
    };
    // The sterol template has no chains for the graph matcher to find, so it
    // is compared whole. `canonical` has been through `canonicalize`, so the
    // reference has to be as well — comparing against the raw emitted form
    // never matched.
    if canonicalize(&generate_smiles("ST").unwrap_or_default()).as_deref() == Some(canonical) {
        return vec!["ST".to_string()];
    }

    let main = connected_components(&mol, &HashSet::new())
        .into_iter()
        .filter(|component| {
            !component
                .iter()
                .any(|&i| mol.atom(AtomIdx(i as u32)).wildcard)
        })
        .max_by_key(Vec::len);
    let Some(main) = main else { return Vec::new() };
    let main_set: HashSet<usize> = main.into_iter().collect();
    let mut out = sphingoid_graph_candidates(&mol, &cx, &main_set);
    let cuts: Vec<(usize, usize)> = mol
        .bonds()
        .filter_map(|(_, bond)| {
            let a = bond.atom1.0 as usize;
            let b = bond.atom2.0 as usize;
            if !main_set.contains(&a) || !main_set.contains(&b) || !single_bond(bond.order) {
                return None;
            }
            match (element(&mol, a), element(&mol, b)) {
                (Element::C, Element::O | Element::N) => Some((a, b)),
                (Element::O | Element::N, Element::C) => Some((b, a)),
                _ => None,
            }
        })
        .collect();

    for template in graph_templates() {
        for chosen in combinations(&cuts, template.chains) {
            if let Some(mut names) = match_graph_template(&mol, &cx, &main_set, template, &chosen) {
                out.append(&mut names);
            }
        }
    }
    out
}

fn sphingoid_graph_candidates(mol: &Molecule, cx: &Cx, main: &HashSet<usize>) -> Vec<String> {
    let mut out = Vec::new();
    for (n_idx, atom) in mol.atoms() {
        if atom.element != Element::N {
            continue;
        }
        let n = n_idx.0 as usize;
        let carbon_neighbors: Vec<usize> = mol
            .neighbors(n_idx)
            .filter_map(|(neighbor, _)| {
                let i = neighbor.0 as usize;
                (main.contains(&i) && element(mol, i) == Element::C).then_some(i)
            })
            .collect();
        for &c2 in &carbon_neighbors {
            let base_neighbors: Vec<usize> = mol
                .neighbors(AtomIdx(c2 as u32))
                .filter_map(|(neighbor, _)| {
                    let i = neighbor.0 as usize;
                    (i != n && element(mol, i) == Element::C).then_some(i)
                })
                .collect();
            if base_neighbors.len() != 2 {
                continue;
            }
            for &c1 in &base_neighbors {
                let Some(&c3) = base_neighbors.iter().find(|&&i| i != c1) else {
                    continue;
                };
                let tail_o = mol.neighbors(AtomIdx(c1 as u32)).find_map(|(neighbor, _)| {
                    let i = neighbor.0 as usize;
                    (i != c2 && element(mol, i) == Element::O).then_some(i)
                });
                let Some(tail_o) = tail_o else { continue };

                let acyl_start = carbon_neighbors.iter().copied().find(|&candidate| {
                    candidate != c2
                        && mol
                            .neighbors(AtomIdx(candidate as u32))
                            .any(|(neighbor, bond)| {
                                element(mol, neighbor.0 as usize) == Element::O
                                    && mol.bond(bond).order == BondOrder::Double
                            })
                });
                let has_n_acyl = acyl_start.is_some();
                let removed = acyl_start
                    .map(|carbonyl| HashSet::from([normalized_edge(n, carbonyl)]))
                    .unwrap_or_default();
                let components = connected_components(mol, &removed);
                let base_component = components.iter().find(|part| part.contains(&c1));
                let Some(base_component) = base_component else {
                    continue;
                };
                let base_allowed: HashSet<usize> = base_component.iter().copied().collect();

                let mut tail_path = Vec::new();
                let mut carbon_allowed = base_allowed.clone();
                carbon_allowed.remove(&c1);
                carbon_allowed.remove(&c2);
                longest_carbon_path(mol, c3, &carbon_allowed, &mut Vec::new(), &mut tail_path);
                if tail_path.is_empty() {
                    continue;
                }
                let mut path = vec![c1, c2];
                path.extend(tail_path);
                let carbon_set: HashSet<usize> = path.iter().copied().collect();
                let tail_atoms = component_without(mol, tail_o, &carbon_set, &base_allowed);
                let Some((tail_mol, _)) = induced_molecule(mol, &tail_atoms) else {
                    continue;
                };
                let tail_signature = canonical_smiles(&tail_mol);

                for (class, tail, expected_acyl) in SPHINGO {
                    if *expected_acyl != has_n_acyl {
                        continue;
                    }
                    let Ok(probe) = parse(tail) else { continue };
                    if canonical_smiles(&probe) != tail_signature {
                        continue;
                    }
                    let Some(mut base) =
                        graph_sphingoid_chain(mol, cx, &base_allowed, &path, n, tail_o)
                    else {
                        continue;
                    };
                    base = match add_floating_mods(mol, cx, base, &carbon_set) {
                        Some(base) => base,
                        None => continue,
                    };
                    let base_token = format_chain(&base);
                    if let Some(acyl_start) = acyl_start {
                        let Some(acyl_part) = components.iter().find(|p| p.contains(&acyl_start))
                        else {
                            continue;
                        };
                        let Some((acyl, carbons)) = graph_chain(mol, cx, acyl_part, acyl_start)
                        else {
                            continue;
                        };
                        let Some(acyl) = add_floating_mods(mol, cx, acyl, &carbons) else {
                            continue;
                        };
                        out.push(format!("{class} {base_token}/{}", format_chain(&acyl)));
                    } else {
                        out.push(format!("{class} {base_token}"));
                    }
                }
            }
        }
    }
    out
}

fn graph_sphingoid_chain(
    mol: &Molecule,
    cx: &Cx,
    allowed: &HashSet<usize>,
    path: &[usize],
    nitrogen: usize,
    tail_oxygen: usize,
) -> Option<Chain> {
    let carbon_set: HashSet<usize> = path.iter().copied().collect();
    let hydroxyl = |carbon: usize| {
        mol.neighbors(AtomIdx(carbon as u32))
            .any(|(neighbor, bond)| {
                let i = neighbor.0 as usize;
                allowed.contains(&i)
                    && element(mol, i) == Element::O
                    && mol.bond(bond).order == BondOrder::Single
                    && mol.degree(neighbor) == 1
            })
    };
    if path.len() < 3 || !hydroxyl(path[2]) {
        return None;
    }
    let triol = path.get(3).is_some_and(|&c| hydroxyl(c));
    let mut chain = Chain {
        prefix: if triol { "t" } else { "d" },
        carbon: path.len() as u32,
        ..Chain::default()
    };
    for (i, pair) in path.windows(2).enumerate() {
        if mol
            .bond_between(AtomIdx(pair[0] as u32), AtomIdx(pair[1] as u32))?
            .1
            .order
            == BondOrder::Double
        {
            chain
                .db
                .push((i as u32 + 1, double_bond_geometry(mol, path, i)));
        }
    }
    let mut handled = HashSet::new();
    for (pos, &carbon) in path.iter().enumerate() {
        for (neighbor, bond_idx) in mol.neighbors(AtomIdx(carbon as u32)) {
            let root = neighbor.0 as usize;
            if carbon_set.contains(&root)
                || !allowed.contains(&root)
                || handled.contains(&root)
                || root == nitrogen
                || root == tail_oxygen
            {
                continue;
            }
            let branch_atoms = component_without(mol, root, &carbon_set, allowed);
            handled.extend(branch_atoms.iter().copied());
            if (pos == 2 || (triol && pos == 3))
                && branch_atoms.len() == 1
                && element(mol, root) == Element::O
            {
                continue;
            }
            let path_attachments = branch_atoms
                .iter()
                .flat_map(|&atom| mol.neighbors(AtomIdx(atom as u32)))
                .filter(|(atom, _)| carbon_set.contains(&(atom.0 as usize)))
                .count();
            if path_attachments == 2 && branch_atoms.len() == 1 && element(mol, root) == Element::O
            {
                let other = mol.neighbors(neighbor).find_map(|(atom, _)| {
                    let atom = atom.0 as usize;
                    (atom != carbon && carbon_set.contains(&atom)).then_some(atom)
                })?;
                let other_pos = path.iter().position(|x| *x == other)?;
                chain.rings.push((
                    pos.min(other_pos) as u32 + 1,
                    pos.max(other_pos) as u32 + 1,
                    true,
                ));
                continue;
            }
            chain.mods.push((
                pos as u32 + 1,
                branch_abbreviation(mol, &branch_atoms, mol.bond(bond_idx).order)?,
            ));
        }
    }
    let (unlocalized, hidden) = cx.scaffold_atoms(&carbon_set);
    chain.unlocalized_db = unlocalized;
    chain.carbon += hidden;
    for _ in 0..unlocalized {
        chain.db.pop();
    }
    chain.db.sort();
    chain.mods.sort();
    Some(chain)
}

fn graph_templates() -> &'static [GraphTemplate] {
    static TEMPLATES: OnceLock<Vec<GraphTemplate>> = OnceLock::new();
    TEMPLATES.get_or_init(|| {
        [
            ("FA", "FA 8:0"),
            ("AMP-FA", "AMP-FA 8:0"),
            ("CE", "CE 8:0"),
            ("CAR", "CAR 8:0"),
            ("NAE", "NAE 8:0"),
            ("MG", "MG 8:0"),
            ("DG", "DG 8:0/9:0"),
            ("TG", "TG 8:0/9:0/10:0"),
            ("LPC", "LPC 8:0"),
            ("PC", "PC 8:0/9:0"),
            ("LPE", "LPE 8:0"),
            ("PE", "PE 8:0/9:0"),
            ("LPS", "LPS 8:0"),
            ("PS", "PS 8:0/9:0"),
            ("LPG", "LPG 8:0"),
            ("PG", "PG 8:0/9:0"),
            ("LPI", "LPI 8:0"),
            ("PI", "PI 8:0/9:0"),
            ("LPA", "LPA 8:0"),
            ("PA", "PA 8:0/9:0"),
            ("CL", "CL 8:0/9:0/10:0/11:0"),
        ]
        .into_iter()
        .filter_map(|(class, reference)| make_graph_template(class, reference))
        .collect()
    })
}

fn make_graph_template(class: &'static str, reference: &'static str) -> Option<GraphTemplate> {
    let structure = crate::forward::generate_structure(reference)?;
    let mol = parse(&structure.smiles).ok()?;
    let mut removed = HashSet::new();
    let mut boundaries = Vec::new();
    for chain in &structure.chains {
        let carbons: HashSet<usize> = chain.carbons.iter().copied().collect();
        let c1 = *chain.carbons.first()?;
        let outside: Vec<(usize, usize)> = mol
            .neighbors(AtomIdx(c1 as u32))
            .filter_map(|(neighbor, bond)| {
                let n = neighbor.0 as usize;
                let order = mol.bond(bond).order;
                (!carbons.contains(&n)
                    && single_bond(order)
                    && matches!(element(&mol, n), Element::O | Element::N))
                .then_some((c1, n))
            })
            .collect();
        let &(carbon, head) = outside.first()?;
        let edge = normalized_edge(carbon, head);
        removed.insert(edge);
        boundaries.push((chain.sn, carbon, head));
    }
    let components = connected_components(&mol, &removed);
    let chain_starts: HashSet<usize> = boundaries.iter().map(|(_, c, _)| *c).collect();
    let head = components
        .into_iter()
        .find(|part| part.iter().all(|atom| !chain_starts.contains(atom)))?;
    let (head_mol, _) = induced_molecule(&mol, &head)?;
    let sn = boundaries.into_iter().map(|(sn, _, _)| sn).collect();
    Some(GraphTemplate {
        class,
        chains: structure.chains.len(),
        sn,
        head_signature: canonical_smiles(&head_mol),
    })
}

fn match_graph_template(
    mol: &Molecule,
    cx: &Cx,
    main: &HashSet<usize>,
    template: &GraphTemplate,
    cuts: &[(usize, usize)],
) -> Option<Vec<String>> {
    let removed: HashSet<(usize, usize)> = cuts
        .iter()
        .map(|&(carbon, head)| normalized_edge(carbon, head))
        .collect();
    let parts: Vec<Vec<usize>> = connected_components(mol, &removed)
        .into_iter()
        .filter(|part| part.iter().any(|atom| main.contains(atom)))
        .collect();
    if parts.len() != template.chains + 1 {
        return None;
    }

    for head in &parts {
        let (head_mol, _) = induced_molecule(mol, head)?;
        if canonical_smiles(&head_mol) != template.head_signature {
            continue;
        }
        let head_atoms: HashSet<usize> = head.iter().copied().collect();
        let mut chain_tokens = Vec::new();
        for &(chain_start, input_head) in cuts {
            if !head_atoms.contains(&input_head) {
                return None;
            }
            let part = parts
                .iter()
                .find(|part| part.contains(&chain_start) && !part.contains(&input_head))?;
            let (chain, carbons) = graph_chain(mol, cx, part, chain_start)?;
            chain_tokens.push(format_chain(&add_floating_mods(mol, cx, chain, &carbons)?));
        }
        let sep = if cx.swappable { "_" } else { "/" };
        let names = permutations(&chain_tokens)
            .into_iter()
            .filter_map(|tokens| {
                let slots = template
                    .sn
                    .iter()
                    .copied()
                    .zip(tokens.into_iter().map(Some))
                    .collect();
                let (_, tokens) = join_slots(slots, sep)?;
                Some(format!("{} {tokens}", template.class))
            })
            .collect();
        return Some(names);
    }
    None
}

fn graph_chain(
    mol: &Molecule,
    cx: &Cx,
    component: &[usize],
    start: usize,
) -> Option<(Chain, HashSet<usize>)> {
    let allowed: HashSet<usize> = component.iter().copied().collect();
    let mut path = Vec::new();
    longest_carbon_path(mol, start, &allowed, &mut Vec::new(), &mut path);
    if path.is_empty() {
        return None;
    }
    let carbon_set: HashSet<usize> = path.iter().copied().collect();
    let mut chain = Chain {
        carbon: path.len() as u32,
        ..Chain::default()
    };

    let c1 = path[0];
    let acyl = mol.neighbors(AtomIdx(c1 as u32)).any(|(neighbor, bond)| {
        element(mol, neighbor.0 as usize) == Element::O && mol.bond(bond).order == BondOrder::Double
    });
    if !acyl {
        chain.prefix = if path.len() > 1
            && mol
                .bond_between(AtomIdx(path[0] as u32), AtomIdx(path[1] as u32))?
                .1
                .order
                == BondOrder::Double
        {
            "P-"
        } else {
            "O-"
        };
    }

    for (i, pair) in path.windows(2).enumerate() {
        let bond = mol
            .bond_between(AtomIdx(pair[0] as u32), AtomIdx(pair[1] as u32))?
            .1;
        if bond.order == BondOrder::Double {
            let geom = double_bond_geometry(mol, &path, i);
            chain.db.push((i as u32 + 1, geom));
        }
    }

    // Extra carbon-carbon edges close Table 1B rings.
    for (_, bond) in mol.bonds() {
        let a = bond.atom1.0 as usize;
        let b = bond.atom2.0 as usize;
        let (Some(pa), Some(pb)) = (
            path.iter().position(|x| *x == a),
            path.iter().position(|x| *x == b),
        ) else {
            continue;
        };
        if pa.abs_diff(pb) > 1 {
            chain
                .rings
                .push((pa.min(pb) as u32 + 1, pa.max(pb) as u32 + 1, false));
        }
    }

    let mut handled = HashSet::new();
    for (pos, &carbon) in path.iter().enumerate() {
        for (neighbor, bond_idx) in mol.neighbors(AtomIdx(carbon as u32)) {
            let root = neighbor.0 as usize;
            if carbon_set.contains(&root) || !allowed.contains(&root) || handled.contains(&root) {
                continue;
            }
            let branch_atoms = component_without(mol, root, &carbon_set, &allowed);
            handled.extend(branch_atoms.iter().copied());
            let path_attachments = branch_atoms
                .iter()
                .flat_map(|&atom| mol.neighbors(AtomIdx(atom as u32)))
                .filter(|(atom, _)| carbon_set.contains(&(atom.0 as usize)))
                .count();
            if path_attachments == 2 && branch_atoms.len() == 1 && element(mol, root) == Element::O
            {
                let other = mol.neighbors(neighbor).find_map(|(atom, _)| {
                    let atom = atom.0 as usize;
                    (atom != carbon && carbon_set.contains(&atom)).then_some(atom)
                })?;
                let other_pos = path.iter().position(|x| *x == other)?;
                chain.rings.push((
                    pos.min(other_pos) as u32 + 1,
                    pos.max(other_pos) as u32 + 1,
                    true,
                ));
                continue;
            }
            if pos == 0
                && element(mol, root) == Element::O
                && mol.bond(bond_idx).order == BondOrder::Double
            {
                continue;
            }
            let abbreviation = branch_abbreviation(mol, &branch_atoms, mol.bond(bond_idx).order)?;
            chain.mods.push((pos as u32 + 1, abbreviation));
        }
    }

    let (unlocalized, hidden) = cx.scaffold_atoms(&carbon_set);
    chain.unlocalized_db = unlocalized;
    chain.carbon += hidden;
    for _ in 0..unlocalized {
        chain.db.pop();
    }
    chain.db.sort();
    chain.mods.sort();
    chain.rings.sort();
    Some((chain, carbon_set))
}

fn add_floating_mods(
    mol: &Molecule,
    cx: &Cx,
    mut chain: Chain,
    carbons: &HashSet<usize>,
) -> Option<Chain> {
    for (floating, sites) in &cx.m {
        if !sites.iter().all(|site| carbons.contains(site)) {
            continue;
        }
        let wildcard = AtomIdx(*floating as u32);
        if !mol.atom(wildcard).wildcard {
            return None;
        }
        let (root, bond) = mol.neighbors(wildcard).next()?;
        let allowed: HashSet<usize> = connected_components(mol, &HashSet::new())
            .into_iter()
            .find(|part| part.contains(floating))?
            .into_iter()
            .filter(|atom| atom != floating)
            .collect();
        let atoms: Vec<usize> = allowed.iter().copied().collect();
        chain
            .mods
            .push((0, branch_abbreviation(mol, &atoms, mol.bond(bond).order)?));
        let _ = root;
    }
    chain.mods.sort();
    Some(chain)
}

impl Cx<'_> {
    fn scaffold_atoms(&self, atoms: &HashSet<usize>) -> (usize, u32) {
        let mut consumed = 0;
        for &(terms, sum) in &self.equations {
            if consumed + terms > self.sg.len() {
                break;
            }
            let markers = &self.sg[consumed..consumed + terms];
            consumed += terms;
            if markers.iter().all(|atom| atoms.contains(atom)) {
                return (terms - 1, (sum - terms) as u32);
            }
        }
        (0, 0)
    }
}

fn branch_abbreviation(
    mol: &Molecule,
    atoms: &[usize],
    attachment: BondOrder,
) -> Option<&'static str> {
    let (branch, _) = induced_molecule(mol, atoms)?;
    let signature = canonical_smiles(&branch);
    SUBSTITUENTS.iter().find_map(|(abbr, text)| {
        let probe = parse(&format!("C({text})")).ok()?;
        let root = 0usize;
        let (neighbor, bond) = probe.neighbors(AtomIdx(root as u32)).next()?;
        let rest: Vec<usize> = (0..probe.atom_count())
            .filter(|atom| *atom != root)
            .collect();
        let (probe_branch, _) = induced_molecule(&probe, &rest)?;
        (same_bond_kind(attachment, probe.bond(bond).order)
            && neighbor.0 > 0
            && canonical_smiles(&probe_branch) == signature)
            .then_some(*abbr)
    })
}

fn longest_carbon_path(
    mol: &Molecule,
    atom: usize,
    allowed: &HashSet<usize>,
    current: &mut Vec<usize>,
    best: &mut Vec<usize>,
) {
    if !allowed.contains(&atom) || element(mol, atom) != Element::C || current.contains(&atom) {
        return;
    }
    current.push(atom);
    if current.len() > best.len() {
        *best = current.clone();
    }
    for (neighbor, _) in mol.neighbors(AtomIdx(atom as u32)) {
        longest_carbon_path(mol, neighbor.0 as usize, allowed, current, best);
    }
    current.pop();
}

fn double_bond_geometry(mol: &Molecule, path: &[usize], bond_pos: usize) -> Option<char> {
    if bond_pos == 0 || bond_pos + 2 >= path.len() {
        return None;
    }
    let before = directed_marker(mol, path[bond_pos - 1], path[bond_pos])?;
    let after = directed_marker(mol, path[bond_pos + 1], path[bond_pos + 2])?;
    Some(if before == after { 'E' } else { 'Z' })
}

fn directed_marker(mol: &Molecule, from: usize, to: usize) -> Option<char> {
    let (_, bond) = mol.bond_between(AtomIdx(from as u32), AtomIdx(to as u32))?;
    match bond.order {
        BondOrder::Up => Some(if bond.atom1.0 as usize == from {
            '/'
        } else {
            '\\'
        }),
        BondOrder::Down => Some(if bond.atom1.0 as usize == from {
            '\\'
        } else {
            '/'
        }),
        _ => None,
    }
}

fn component_without(
    mol: &Molecule,
    start: usize,
    blocked: &HashSet<usize>,
    allowed: &HashSet<usize>,
) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut stack = vec![start];
    while let Some(atom) = stack.pop() {
        if blocked.contains(&atom) || !allowed.contains(&atom) || !seen.insert(atom) {
            continue;
        }
        stack.extend(
            mol.neighbors(AtomIdx(atom as u32))
                .map(|(neighbor, _)| neighbor.0 as usize),
        );
    }
    seen.into_iter().collect()
}

fn induced_molecule(mol: &Molecule, atoms: &[usize]) -> Option<(Molecule, HashMap<usize, usize>)> {
    let mut atoms = atoms.to_vec();
    atoms.sort_unstable();
    let set: HashSet<usize> = atoms.iter().copied().collect();
    let mut builder = MoleculeBuilder::new();
    let mut mapping = HashMap::new();
    for old in atoms {
        let new = builder.add_atom(mol.atom(AtomIdx(old as u32)).clone());
        mapping.insert(old, new.0 as usize);
    }
    for (old_bond, bond) in mol.bonds() {
        let a = bond.atom1.0 as usize;
        let b = bond.atom2.0 as usize;
        if set.contains(&a) && set.contains(&b) {
            let new_bond = builder
                .add_bond(
                    AtomIdx(*mapping.get(&a)? as u32),
                    AtomIdx(*mapping.get(&b)? as u32),
                    bond.order,
                )
                .ok()?;
            if let Some(direction) = mol.bond_direction(old_bond) {
                builder.set_bond_direction(new_bond, direction);
            }
        }
    }
    Some((builder.build(), mapping))
}

/// Every atom reachable from `start` without crossing a `removed` bond.
pub(crate) fn reachable_from(
    mol: &Molecule,
    start: usize,
    removed: &HashSet<(usize, usize)>,
) -> HashSet<usize> {
    let mut seen = HashSet::new();
    let mut stack = vec![start];
    while let Some(atom) = stack.pop() {
        if !seen.insert(atom) {
            continue;
        }
        for (neighbor, _) in mol.neighbors(AtomIdx(atom as u32)) {
            let neighbor = neighbor.0 as usize;
            if !removed.contains(&normalized_edge(atom, neighbor)) {
                stack.push(neighbor);
            }
        }
    }
    seen
}

fn connected_components(mol: &Molecule, removed: &HashSet<(usize, usize)>) -> Vec<Vec<usize>> {
    let mut seen: HashSet<usize> = HashSet::new();
    let mut out = Vec::new();
    for start in 0..mol.atom_count() {
        if seen.contains(&start) {
            continue;
        }
        let mut part: Vec<usize> = reachable_from(mol, start, removed).into_iter().collect();
        part.sort_unstable();
        seen.extend(&part);
        out.push(part);
    }
    out
}

fn combinations<T: Copy>(items: &[T], count: usize) -> Vec<Vec<T>> {
    fn visit<T: Copy>(
        items: &[T],
        count: usize,
        at: usize,
        cur: &mut Vec<T>,
        out: &mut Vec<Vec<T>>,
    ) {
        if cur.len() == count {
            out.push(cur.clone());
            return;
        }
        for i in at..items.len() {
            cur.push(items[i]);
            visit(items, count, i + 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    visit(items, count, 0, &mut Vec::new(), &mut out);
    out
}

fn permutations<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
    fn visit<T: Clone>(items: &mut [T], at: usize, out: &mut Vec<Vec<T>>) {
        if at == items.len() {
            out.push(items.to_vec());
            return;
        }
        for i in at..items.len() {
            items.swap(at, i);
            visit(items, at + 1, out);
            items.swap(at, i);
        }
    }
    let mut items = items.to_vec();
    let mut out = Vec::new();
    visit(&mut items, 0, &mut out);
    out
}

fn normalized_edge(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

fn element(mol: &Molecule, atom: usize) -> Element {
    mol.atom(AtomIdx(atom as u32)).element
}

pub(crate) fn single_bond(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
}

fn same_bond_kind(a: BondOrder, b: BondOrder) -> bool {
    (single_bond(a) && single_bond(b)) || a == b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round(name: &str) -> Option<String> {
        parse_smiles(&generate_smiles(name)?)
    }

    #[test]
    fn plain_chains_come_back_exactly() {
        for name in [
            "FA 18:0",
            "FA 18:1(9)",
            "FA 18:1(9Z)",
            "FA 18:2(9Z,12Z)",
            "PC 16:0/18:1(9Z)",
            "TG 16:0/18:1(9)/18:2(9,12)",
            "CL 16:0/18:1(9)/16:0/18:1(9)",
            "CE 18:1(9Z)",
            "CAR 18:1(9)",
            "NAE 20:4(5,8,11,14)",
            "AMP-FA 20:4(5,8,11,14);15OH",
            "MG 18:1(9)",
            "LPC 18:1(9)",
            "PC O-16:0/18:1(9)",
            "PC P-16:0/18:1(9)",
            "ST",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    /// Undetermined structural features remain undetermined after a round trip.
    #[test]
    fn undetermined_things_stay_undetermined() {
        for name in [
            "FA 18:1",               // unlocalized double bond -> Sg:
            "FA 20:4",               // four of them
            "FA 18:2(9)",            // one placed, one not
            "FA 18:0;OH",            // unlocalized modification -> m:
            "FA 18:1;OH",            // both at once
            "FA 20:3(5,8,11);(OH)2", // two unlocalized groups
            "PC 16:0_18:1(9)",       // unresolved sn assignment
            "TG 16:0_18:1(9)_18:0",  // three-way
            "CL 16:0_18:1(9)_16:0_18:1(9)",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    #[test]
    fn sphingoid_bases_recover_their_prefix_and_n_acyl() {
        for name in [
            "Cer d18:1(4)/16:0",
            "SM d18:1(4)/16:0",
            "HexCer d18:1(4)/16:0",
            "Sph d18:1(4)",
            "S1P d18:1(4)",
            "Cer d18:1(4)/16:0;9Ep",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }
    }

    #[test]
    fn table_1a_groups_and_table_1b_rings_survive() {
        for name in [
            "FA 18:0;5Me",
            "FA 16:0;3Me,7Me,11Me,15Me",
            "FA 18:0;5Br",
            "FA 18:0;5CN",
            "FA 18:1(9);12NO2",
            "FA 18:1(9);(NO2)",
            "FA 20:3(5Z,13E);11OH,15OH;9oxo",
            "FA 18:0;9Ep",
            "FA 19:0;[11-13cy3:0]",
            "FA 19:0;[9-11cy3:1(9)]",
        ] {
            assert_eq!(round(name).as_deref(), Some(name), "{name}");
        }

        // A cyclopentene ring is symmetric about the carbon it hangs from, so
        // a double bond at 15-16 and one at 17-18 are the same molecule. The
        // reverse direction returns one spelling for both; it regenerates the
        // same structure, which is the property that matters.
        let gorlic = "FA 18:1(6Z);[14-18cy5:1(15)]";
        let recovered = round(gorlic).expect("should resolve");
        assert_eq!(
            generate_smiles(&recovered).and_then(|s| canonicalize(&s)),
            generate_smiles(gorlic).and_then(|s| canonicalize(&s)),
            "{gorlic} came back as {recovered}, a different structure"
        );
    }

    /// Not a general SMILES parser, and honest about it.
    #[test]
    fn foreign_structures_are_refused() {
        for smi in [
            "",
            "c1ccccc1",                                 // benzene
            "OC(=O)CCCCCCCCCCCCCCCCC.O",                // a stearate with something extra
            "C(COP(=O)([O-])OCC[N+](C)(C)C)(OC(=O)CCC", // truncated
            "OC(=O)CCCC[Se]CCCC",                       // element this crate never writes
        ] {
            assert_eq!(parse_smiles(smi), None, "{smi:?} should be refused");
        }
    }

    /// A returned name must regenerate the same canonical structure it read.
    /// Tampering with a structure therefore has only two honest outcomes: a
    /// name for the *tampered* structure, or nothing at all — never the name
    /// of the structure it was derived from.
    #[test]
    fn a_tampered_structure_never_returns_the_original_name() {
        let name = "PC 16:0/18:1(9)";
        let original = generate_smiles(name).unwrap();

        for tampered in [
            original.replacen("CCCCCCCC=C", "CCCCCCC=C", 1), // move the double bond
            original.replacen("COC(=O)CCCC", "COC(=O)CCC", 1), // shorten a chain
            original.replace("[N+](C)(C)C", "[N+](C)(C)CC"), // alter the headgroup
        ] {
            assert_ne!(tampered, original, "tamper should have changed something");
            match parse_smiles(&tampered) {
                None => {}
                Some(recovered) => {
                    assert_ne!(recovered, name, "{tampered}: returned the untampered name");
                    let regenerated = generate_smiles(&recovered)
                        .and_then(|s| canonicalize(&s))
                        .expect("recovered name should regenerate");
                    assert_eq!(
                        regenerated,
                        canonicalize(&tampered).unwrap(),
                        "{tampered}: recovered name is not equivalent"
                    );
                }
            }
        }
    }

    /// A chain read back and written out again must be a fixed point: if the
    /// first trip canonicalized anything, the second must change nothing.
    #[test]
    fn the_round_trip_is_idempotent() {
        for name in [
            "FA 18:0;ep(5)",
            "FA 18:0;OH(5)",
            "PC 16:0_18:1",
            "DG 16:0/0:0",
        ] {
            let once = round(name).unwrap_or_else(|| panic!("{name} should resolve"));
            let twice = round(&once).unwrap_or_else(|| panic!("{once} should resolve"));
            assert_eq!(once, twice, "{name} did not settle after one trip");
        }
    }
}
